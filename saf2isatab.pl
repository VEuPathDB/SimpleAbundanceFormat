#!/usr/bin/env perl
# -*- mode: cperl -*- 


#
# usage: ./saf2isatab.pl config.yaml saf-data.tsv [ saf-data2.tsv ... ]
#
# options:
#
#   --entities entities.yaml  (tell it what entities we are using)
#   --skip-defaults           (don't load the default SAF column configurations)
#   --novalidate              (don't check that ontology term IDs are in the owl files - see next)
#   --owl_file_url URL --owl_file_url URL2 ...
#                             (one or more URLs to raw OWL file downloads, e.g. github links)
#
#
# note:
# * can handle tsv or csv data
# * if using multiple data files
#   - all files must contain all the mandatory columns
#   - all files must contain at least one shared root entity
#

#
# TO DO:
#
#- make column 'requiredness' entity-dependent somehow - so when using multiple files, there
#  are required columns IF you have the entity_ID column already defined, so maybe it becomes
#  { true, false, entity_dependent }?
#  - how to make that work with __AUTO__ entity IDs? - maybe we don't and keep using required=false for now
#  - make_auto_entity_id should return '' if all column values are ''
#
#- use ordered hashrefs to preserve order
#
#- split and rejoin non-term multivalued values to change delimiter if needed
#
#- do range checking on lat/long and date format checks
#
#- check for "overwrite" errors: e.g. a collection can occur on multiple rows - its data ought to be identical
#  currently we only load the first row and don't check subsequent rows.
#
#- collect and report all errors rather than die for each one
#
#- check study_protocols for duplicates after appending defaults
#
#- specify format for comments (e.g. semicolon delimited, colon prefixed?) and load them separately according to prefix
#



use strict;
use warnings;

use Getopt::Long;
use YAML::XS qw/LoadFile/;
use Text::CSV::Hashify; # previously we appended a version number 0.11
use Bio::Parser::ISATab;
use FindBin;
use Hash::Merge::Simple qw/merge/;
# use Data::Dumper;
use Storable qw/dclone/;
use JSON;
use Scalar::Util qw/looks_like_number/;
use Tie::Hash::Indexed;
use LWP::Simple;

#
# NOTE ABOUT YAML: LoadFile will autodetect booleans, so 'true'/'false' values
# are loaded as Perl booleans - although they don't really exist, so just 1 and undefined.
#
my $default_config_filename = $FindBin::Bin."/default-config.yaml";
my $defaultConfig = LoadFile($default_config_filename);

my $default_input_delimiter = ';';
my $default_isatab_delimiter = ';';
my $output_dir = 'temp-isatab';
my $entities_filename = $FindBin::Bin."/entities.yaml";
my $skip_defaults;

my @owl_file_urls =
  (
   'https://github.com/VEuPathDB/ApiCommonData/raw/master/Load/ontology/popbio/popbio.owl',
   'https://github.com/VEuPathDB/ApiCommonData/raw/master/Load/ontology/popbio/popbio_taxonomy.owl',
  );
my $validate = 1;

GetOptions("output_directory|output-directory=s" => \$output_dir,
	   "entities=s" => \$entities_filename,
	   "skip-defaults|skip_defaults" => \$skip_defaults,
	   "owl_file_url=s" => \@owl_file_urls,
	   "validate!" => \$validate,
	  );

my ($config_filename, @saf_filenames) = @ARGV;

die "must provide two or more filenames on commandline: config_file saf_data_file(s)\n"
    unless ($config_filename && -e $config_filename && @saf_filenames);

my $valid_term_ids = $validate ? get_valid_term_ids(@owl_file_urls) : {};

my @DEFERRED_ERRORS;

my $userConfig = LoadFile($config_filename);
die "problem reading '$config_filename'\n" unless $userConfig;

# merge the user config into the default config
my $config = $skip_defaults ? $userConfig : merge $defaultConfig, $userConfig;

# pull out the column config
my $column_config = $config->{columns};
# print Dumper($column_config); exit;

# a subset of the config is also the basis of the ISA-Tab study datastructure
# the ISA-Tab writer ignores the extra config information
my $study = $config;

# the YAML merge can't handle arrays, so we append the default_study_protocols onto the main protocols
# (with no checks for duplicates...)
push @{$study->{study_protocols}}, @{$config->{default_study_protocols}} if ($config->{default_study_protocols});

# load the entity graph
my $entities = LoadFile($entities_filename);

die "FATAL ERROR: there must be exactly one root entity (in '$entities_filename')\n"
  if (@$entities != 1);
die "FATAL ERROR: root entity (in '$entities_filename') must be a material type entity\n"
  unless ($entities->[0]{type} eq 'material');

my $root_entity = $entities->[0];
my $flat_entities = flatten($root_entity);
my $entity_to_parent = {};
make_parent_lookup($root_entity, $entity_to_parent);

# make sure column config contains required info
validate_config($config, $flat_entities);
check_config_for_placeholder_strings($config);

if ($validate) {
  # perform deep/recursive validation of term IDs keyed by '*term_accession_number'
  validate_term_accession_numbers($entities, $entities_filename);
  validate_term_accession_numbers($config, $config_filename);
}

foreach my $saf_filename (@saf_filenames) {
  die "FATAL ERROR: data file '$saf_filename' not found\n" unless (-e $saf_filename);

  # load the actual data
  my $hashified = Text::CSV::Hashify->new({
					   file => $saf_filename,
					   key => 'sample_ID',
					   format => 'aoh',
					   sep_char => $saf_filename =~ /\.csv$/ ? ',' : "\t",
					  });

  # an arrayref of the column headings in the input file
  my $column_keys = $hashified->fields;

  # all the rows
  my $rows = $hashified->all;

  # make sure no required columns are missing from column_keys
  # and warn about any unconfigured columns
  $column_keys = validate_columns($column_keys, $column_config);

  # append mandatory columns that have default values
  $column_keys = add_mandatory_columns($column_keys, $column_config);

  $study->{study_file_name} = 's_samples.txt';

  my $study_assays = $study->{study_assays} = [];
  my $study_protocols = $study->{study_protocols} //= [];

  foreach my $row (@$rows) {
    # add material entities (descending the entity graph into assay entities also)
    check_row_for_placeholder_strings($row, $config);
    add_material($root_entity, $row, $study, $column_keys, $config);
  }

}

if (@DEFERRED_ERRORS) {
  die "Cannot proceed with writing ISA-Tab due to the following data errors:\n\n @DEFERRED_ERRORS\n";
}
#
# write the ISA-Tab!
#
warn "Writing ISA-Tab to '$output_dir'...\n";
my $isa_writer = Bio::Parser::ISATab->new(directory => $output_dir, protocols_first=>1);
$isa_writer->write( { ontologies => [], studies => [ $study ] } );



sub add_material {
  my ($entity, $row, $isaref, $column_keys, $config_and_study) = @_;

  my $column_config = $config_and_study->{columns};
  # figure out an ID for this entity
  # this is simple if all *_ID fields are mandatory
  # otherwise we'll need to generate default IDs (which we can add later if providing, e.g. location_ID is a hassle)

  my $id_column_name = find_id_column_name($column_config, $entity);

  die "FATAL ERROR: couldn't find ID column name for $entity->{name}\n"
    unless ($id_column_name);

  my $entity_id = $row->{$id_column_name};
  # warn "id column is $id_column_name and got $entity_id\n";

  if (!defined $entity_id &&
      defined $column_config->{$id_column_name}{default} &&
      $column_config->{$id_column_name}{default} eq '__AUTO__') {
    $entity_id = make_auto_entity_id($entity, $row, $column_keys, $config_and_study, $entity_to_parent);
    # store it in the row data (nasty side effect!)
    # so that make_auto_entity_id can use auto-generated parent IDs
    $row->{$id_column_name} = $entity_id;
  }

  # if there's no entity_id, silently skip processing this entity
  # (typical use case: leave sample_ID column empty for collections that found nothing at all)
  return unless (defined $entity_id && $entity_id ne '');

  # make a hashref to put the new material nodes in, e.g.
  # $study->{sources}{some_source_id} = ...
  $isaref->{$entity->{isa_key}} //= ordered_hashref();
  my $new_isaref = $isaref->{$entity->{isa_key}}{$entity_id};

  # only make a node once (i.e. don't reprocess a location multiple times for each collection and sample)
  if (!defined $new_isaref) {
    $new_isaref = $isaref->{$entity->{isa_key}}{$entity_id} //= { material_type => $entity->{isa_material_type} };

    add_column_data($new_isaref, $column_keys, $row, $entity, $config_and_study);
  }

  # recurse down entity tree
  foreach my $child_entity (@{$entity->{children}}) {
    # check material or assay
    if ($child_entity->{type} eq 'material') {
      add_material($child_entity, $row, $new_isaref, $column_keys, $config_and_study);
    } elsif ($child_entity->{type} eq 'assay') {
      add_assay($child_entity, $row, $new_isaref, $column_keys, $config_and_study);
    }
  }

}


#
# add data in $column_keys from $row to hang off $isaref (in the nascent ISA-Tab data structure)
# e.g. characteristics, comments, protocol ref
#
sub add_column_data {
  my ($isaref, $column_keys, $row, $entity, $config_and_study) = @_;
  my $column_config = $config_and_study->{columns};

  my $relevant_columns = [ grep { $column_config->{$_}{describes} eq $entity->{name} } @$column_keys ];

  foreach my $column (@$relevant_columns) {
    my $col_config = $column_config->{$column};
    my $col_term = $col_config->{column_term};
    my $value = $row->{$column};

    if (!defined $value || $value eq '') {
      $value = $col_config->{default} // '';
    }

    next unless (defined $value && $value ne '');

    # handle the characteristics/variables (not comments or protocols)
    if ($col_term) {
      my $characteristics = ($isaref->{characteristics} //= ordered_hashref())->{"$column (REF:$col_term)"} //= ordered_hashref();

      # if it's a plain text/number/date value then it's a simple case
      # multivalued values can be left as they are
      if ($col_config->{value_type} =~ /^(number|string|date|latitude|longitude)$/) {
	$characteristics->{value} = $value;
	# check if it was allowed
	if ($col_config->{allowed_values}) {
	  push @DEFERRED_ERRORS, "value '$value' not in 'allowed_values' for column '$column'\n"
	    unless (grep { $value eq $_ } @{$col_config->{allowed_values}});
	  # not very efficient with the grep but 'allowed_values' will typically be small
	}

	# check that number values really are numbers, for the appropriate column types
	if ($col_config->{value_type} =~ /^(number|latitude|longitude)$/) {
	  push @DEFERRED_ERRORS, "value '$value' is not numeric in number or lat/long column '$column'\n"
	    unless (looks_like_number($value));
	}

	# now handle units, if provided
	if ($col_config->{unit}) {
	  $characteristics->{unit}{value} = $col_config->{unit};
	  if ($col_config->{unit_term}) {
	    $characteristics->{unit}{term_source_ref} = 'TERM';
	    $characteristics->{unit}{term_accession_number} = $col_config->{unit_term};
	  }
	}

	# handle scale_factor if provided - will always be numbers due to config validation
	if (defined $col_config->{scale_factor}) {
	  $characteristics->{value} *= $col_config->{scale_factor};
	}

	# ontology term values require a lookup from text to term:
      } elsif ($col_config->{value_type} eq 'term') {
	# get the lookup hash (already validated - no need to check success)
	my $lookup = $config->{$col_config->{term_lookup} // 'study_terms'};

	# unfortunately, multivalued term columns make things complicated
	my @values = $col_config->{multivalued} ? split /\s*$col_config->{delimiter}\s*/, $value : ($value);
	my @term_source_refs;
	my @term_accession_numbers;

	foreach my $value (@values) {
	  next unless (defined $value && $value ne ''); # don't do lookups on nothing
	  my $value_term_id = $lookup->{$value};
	  if ($value_term_id) {
	    push @term_source_refs, source_ref($value_term_id);
	    push @term_accession_numbers, $value_term_id;

	    if ($validate && !$valid_term_ids->{$value_term_id}) {
	      push @DEFERRED_ERRORS, sprintf "unvalidated term ID '$value_term_id'; not found in OWL file(s)\n";
	    }
	  } else {
	    push @DEFERRED_ERRORS, sprintf "value '%s' not found in '%s' term lookup\n", $value, $col_config->{term_lookup} // 'study_terms';
	  }
	}
	$characteristics->{value} = join $default_isatab_delimiter, @values;
	$characteristics->{term_source_ref} = join $default_isatab_delimiter, @term_source_refs;
	$characteristics->{term_accession_number} = join $default_isatab_delimiter, @term_accession_numbers;
      }
    } elsif ($col_config->{value_type} eq 'protocol_ref') {
      # always assume splitting as multivalued
      foreach my $p_ref (split /\s*$col_config->{delimiter}\s*/, $value) {
	# check that the protocol ref is in the study_protocols
	my @ok = grep { $_->{study_protocol_name} eq $p_ref } @{$config_and_study->{study_protocols}};
	if (@ok) {
	  $isaref->{protocols}{$p_ref} = ordered_hashref();
	} else {
	  push @DEFERRED_ERRORS, "protocol ref '$p_ref' not found in study_protocols\n";
	}
      }
    } elsif ($col_config->{value_type} eq 'comment') {
      $isaref->{comments}{$column} = $value;
    }
  }
}


sub add_assay {
  my ($entity, $row, $isaref, $column_keys, $config_and_study) = @_;

  # There can be multiple assays for each type, e.g. 'insecticide resistance assay',
  # when there are different protocols used for different columns with the `protocol` property.
  # So we have to process all the columns and figure out which assay they correspond to.

  my $study_assay_measurement_type = $entity->{name};

  # there are two ways to find out what protocol is needed
  # 1. A column such as `species_identification_method` provides the protocol on a row by row basis.
  #    For this we search the column config and make $protocol_column - if there is just one (e.g. `species_identification_method`)
  #    then we'll use that.
  # 2. The assay type (e.g. `insecticide resistance assay`) doesn't have a special *_method column but
  #    individual columns can have a `protocol` attribute, and we use that on a column-wise basis.
  #
  my $protocol_columns = [ grep {
    $column_config->{$_}{describes} eq $entity->{name} &&
      $column_config->{$_}{value_type} eq 'protocol_ref'
    } @$column_keys ];

  my $protocol_column;
  if (@$protocol_columns == 1) {
    $protocol_column = $protocol_columns->[0];
  } elsif (@$protocol_columns > 1) {
    die "FATAL ERROR: didn't expect more than one protocol_ref column for $study_assay_measurement_type: ".join(', ', @$protocol_columns)."\n";
  }

  my $data_columns = [ grep {
    $column_config->{$_}{describes} eq $entity->{name} &&
      $column_config->{$_}{column_term}
  } @$column_keys ];

  foreach my $column (@$data_columns) {
    my $col_config = $column_config->{$column};
    my $col_term = $col_config->{column_term};
    my $value = $row->{$column};
    my $row_protocol_ref = $protocol_column ? $row->{$protocol_column} || $column_config->{$protocol_column}{default} : undef;
    my $column_protocol_ref = $col_config->{protocol};

    die "FATAL ERROR: no row- or column-wise protocol for sample '$row->{sample_ID}' and '$study_assay_measurement_type'\n"
      unless ($row_protocol_ref || $column_protocol_ref);

    # now we need to figure out which study_assay to load the data into
    my $assay_filename = make_assay_filename($study_assay_measurement_type, $row_protocol_ref, $column_protocol_ref);
    my $study_assay = find_or_create_study_assay($config_and_study, $study_assay_measurement_type, $assay_filename);

    # make an assay ID - which will be '' if there is no data for the assay
    my $assay_id = make_auto_entity_id($entity, $row, $column_keys, $config_and_study, $entity_to_parent);

    if ($assay_id) {
      # now add the link sample and assay
      my $assay = ($study_assay->{samples}{$row->{sample_ID}}{assays} //= ordered_hashref())->{$assay_id} //= ordered_hashref();

      # add the column-wise protocol ref explicitly
      $assay->{protocols}{$column_protocol_ref} = ordered_hashref()
	if ($column_protocol_ref);
      # otherwise row-wise protocol ref will be added from the input $row by add_column_data
      add_column_data($assay, $column_keys, $row, $entity, $config_and_study);
    }
  }
}

sub validate_config {
  my ($config, $flat_entities) = @_;

  my $column_config = $config->{columns};
  # check that every column definition has `describes` and `value_type`
  validate_cols($column_config, sub {
		  my $col = shift;
		  (!$col->{describes} || !$col->{value_type});
		}, "column configuration is missing required attributes 'describes' and/or 'value_type' for columns");

  # now check that all the columns `describe` an existing entity
  my %entity_names = map { ($_->{name} => 1) } @$flat_entities;
  validate_cols($column_config, sub {
		  my $col = shift;
		  !$entity_names{$col->{describes}};
		}, "these columns 'describe' entities that do not exist");

  # check that value_type is supported
  my $value_type_validation = qr/^(term|number|string|date|latitude|longitude|id|comment|protocol_ref)$/;
  validate_cols($column_config, sub {
		  my $col = shift;
		  $col->{value_type} !~ $value_type_validation;
		}, "these columns have an unsupported 'value_type'");

  # check that every non-deprecated term|number|string|date|latitude|longitude column has a column_term
  my $value_types_requiring_column_term = qr/^(term|number|string|date|latitude|longitude)$/;
  validate_cols($column_config, sub {
		  my $col = shift;
		  $col->{value_type} =~ $value_types_requiring_column_term &&
		    !$col->{deprecated} &&
		    !$col->{column_term};
		}, "these columns don't have a column_term (e.g. a variable IRI)");

  # validate column_term IDs using OWL file(s)
  if ($validate) {
    my @bad_cols = validate_cols($column_config, sub {
				   my $col = shift;
				   $col->{value_type} =~ $value_types_requiring_column_term &&
				     !$col->{deprecated} &&
				     !$col->{ignore} &&
				     $col->{column_term} &&
				     !$valid_term_ids->{$col->{column_term}};
				 });
    if (@bad_cols > 0) {
      die "FATAL ERROR: the following columns have 'column_term' IDs that were not validated via the OWL file(s):\n".join("\n", map { " $_ (bad ID: $column_config->{$_}{column_term})" } @bad_cols)."\n";
    }
  }

  # check that term columns do not have allowed_values
  validate_cols($column_config, sub {
		  my $col = shift;
		  $col->{value_type} eq 'term' &&
		    defined $col->{allowed_values};
		}, "these 'value_type: term' columns cannot have 'allowed_values'");

  # check that allowed_values are only provided for string and number fields
  # (pretty much a repeat of the above, but we'll keep it for now)
  my $value_types_suitable_for_allowed_values = qr/^(string|number)$/;
  validate_cols($column_config, sub {
		  my $col = shift;
		  defined $col->{allowed_values} &&
		    $col->{value_type} !~ $value_types_suitable_for_allowed_values;
		}, "these columns are not suitable for 'allowed_values'");

  # check that a default value, if given, is in allowed_values, if given
  validate_cols($column_config, sub {
		  my $col = shift;
		  defined $col->{allowed_values} &&
		    $col->{default} &&
		    0 == grep { $_ eq $col->{default} } @{$col->{allowed_values}};
		}, "'default' values for these columns are not in their 'allowed_values'");

  # check that any term_lookup values exist in the $config hash as first level keys
  validate_cols($column_config, sub {
		  my $col = shift;
		  $col->{term_lookup} &&
		    !exists $config->{$col->{term_lookup}};
		}, "these columns have term_lookup values that are not defined in the config file");

  # check that any `protocol` values for assay variables are in the study_protocols
  validate_cols($column_config, sub {
		  my $col = shift;
		  $col->{protocol} &&
		    !grep { $_->{study_protocol_name} eq $col->{protocol} } @{$config->{study_protocols}};
		}, "these columns contain `protocol` values that are not in study_protocols");

  # check that any column with `unit` annotation also has `value_type: number`
  validate_cols($column_config, sub {
		  my $col = shift;
		  $col->{unit} &&
		    $col->{value_type} ne 'number';
		}, "these columns have units but are not number columns");

  # check that scale_factor annotations only apply to number variables
  validate_cols($column_config, sub {
		  my $col = shift;
		  defined $col->{scale_factor} && # scale_factor could be zero in extreme edge cases
		    $col->{value_type} ne 'number';
		}, "these columns have scale_factor but are not number columns");

  # check that scale_factor is a number
  validate_cols($column_config, sub {
		  my $col = shift;
		  defined $col->{scale_factor} && !looks_like_number($col->{scale_factor});
		}, "the scale_factor annotation must be a number for these columns");

  ### the following have side effects!

  # add the default `required: true` to any column that doesn't have it
  map { $_->{required} //= 1 } values %$column_config;

  # add default delimiter for multivalued variables
  map { $_->{delimiter} //= $default_input_delimiter } grep { $_->{multivalued} || (defined $_->{value_type} && $_->{value_type} eq 'protocol_ref') } values %$column_config;
}

# not to be confused with validate_columns below!
# this is just a helper for validate_config
# note that all tests include the `ignore` check
sub validate_cols {
  my ($column_config, $condition, $error_message) = @_;
  my @bad = grep {
    my $col = $column_config->{$_};
    !$col->{ignore} && $condition->($col);
  } keys %$column_config;
  if ($error_message) {
    die "FATAL ERROR: $error_message: " . join(', ', @bad) . "\n" if @bad;
  }
  return @bad;
}

sub validate_term_accession_numbers {
  my ($data, $filename) = @_;
  my @accessions = grep_keys_get_values($data, qr/term_accession_number$/);
  my @bad_accessions = grep { !$valid_term_ids->{$_} } @accessions;
  if (@bad_accessions) {
    die "FATAL ERROR: term_accession_number(s) in '$filename' could not be validated via OWL file(s): ".
      join(", ", @bad_accessions)."\n";
  }
}


sub check_config_for_placeholder_strings {
  my ($config) = @_;
  my $copy = dclone($config);
  my $placeholder_strings = delete $copy->{placeholder_strings};
  my $config_json = encode_json($copy);
  foreach my $placeholder_string (@$placeholder_strings) {
    if ($config_json =~ /$placeholder_string/) {
      push @DEFERRED_ERRORS, "placeholder string '$placeholder_string' was found in your config file\n";
    }
  }
}

sub check_row_for_placeholder_strings {
  my ($row, $config) = @_;
  my $placeholder_strings = $config->{placeholder_strings};
  my $row_json = encode_json($row);
  foreach my $placeholder_string (@$placeholder_strings) {
    if ($row_json =~ /$placeholder_string/) {
      my $row_IDs_only = { map { ($_ => $row->{$_}) } grep { /ID$/ } keys %$row };
      my $row_IDs_json = encode_json($row_IDs_only);
      push @DEFERRED_ERRORS, "placeholder string '$placeholder_string' was found in a row of your data file ($row_IDs_json)\n";
    }
  }
}

sub validate_columns {
  my ($column_keys, $column_config) = @_;

  # make a hash look-up
  my %column_keys = map { ($_ => 1) } @$column_keys;

  # are all the required columns in the input file's header?
  #
  my @missing = grep {
    !$column_config->{$_}{ignore} &&
    $column_config->{$_}{required} &&
    !defined $column_config->{$_}{default} &&
    !$column_keys{$_}
  } keys %$column_config;

  die "FATAL ERROR: the following required columns are missing from the input file: ".join(', ', @missing)."\n"
    if (@missing);

  my @unconfigured = grep {
    !exists $column_config->{$_}
  } @$column_keys;

  warn "WARNING: the following data file columns are not configured (or explicitly ignored) and will be skipped: ".join(', ', @unconfigured)."\n"
    if (@unconfigured);

  # return only the configured and non-ignored columns
  return [ grep { exists $column_config->{$_} && !$column_config->{$_}{ignore} } @$column_keys ];
}

sub add_mandatory_columns {
  my ($column_keys, $column_config) = @_;
  # make a hash look-up
  my %column_keys = map { ($_ => 1) } @$column_keys;
  # work with a copy
  my $result = [ @$column_keys ];

  push @$result,
    grep {
      !$column_keys{$_} && # if we don't already have it
      !$column_config->{$_}{ignore} && # and if we're not ignoring it
      $column_config->{$_}{required} && # and it's required
      defined $column_config->{$_}{default} # and there's a default provided
    } keys %$column_config;

  return $result;
}

sub flatten {
  my ($entity) = @_;
  return [ $entity, map { @{ flatten($_) } } @{$entity->{children}} ];
}


#
# map each node to its immediate parent, in the $lookup hash
#
sub make_parent_lookup {
  my ($entity, $lookup, $parent) = @_;
  # Base case: If the entity is undefined, return.
  return unless defined $entity;

  # If a parent exists, map the current entity to its parent.
  $lookup->{$entity} = $parent if defined $parent;

  # Recursively call this function for each child, passing the current entity as the new parent.
  foreach my $child (@{$entity->{children}}) {
    make_parent_lookup($child, $lookup, $entity);
  }
}

sub make_assay_filename {
  my ($study_assay_measurement_type, $row_protocol_ref, $column_protocol_ref) = @_;

  my $result;
  if ($row_protocol_ref) {
    $result = "a_$study_assay_measurement_type.txt";
  } else {
    $result = "a_$study_assay_measurement_type $column_protocol_ref.txt";
  }
  $result =~ s/\s/_/g;
  return $result;
}


#
# look up in study_assays based on the assay filename
#
sub find_or_create_study_assay {
  my ($config_and_study, $study_assay_measurement_type, $assay_filename) = @_;

  my @existing = grep {
    $_->{study_assay_file_name} eq $assay_filename
  } @{$config_and_study->{study_assays}};

  if (@existing == 1) {
    return $existing[0];
  } elsif (@existing > 1) {
    die "FATAL ERROR: more than one study_assay with same filename\n";
  }

  # now create a new one if necessary
  my $study_assay =
    {
     study_assay_file_name => $assay_filename,
     study_assay_measurement_type => $study_assay_measurement_type,
     study_assay_measurement_type_term_source_ref => source_ref($config_and_study->{study_assay_measurement_type_term_lookup}{$study_assay_measurement_type}),
     study_assay_measurement_type_term_accession_number => $config_and_study->{study_assay_measurement_type_term_lookup}{$study_assay_measurement_type},
     samples => ordered_hashref(),
    };

  push @{$config_and_study->{study_assays}}, $study_assay;
  return $study_assay;
}


#
# generate an entity_id based on the unique signature of all its column values
#
# but return empty string (no ID) if all its columns are empty
#

# we use the size of this hash to generate the entity ID serial number
my %seen_signatures;

sub make_auto_entity_id {
  my ($entity, $row, $column_keys, $config_and_study, $entity_to_parent) = @_;

  my $column_config = $config_and_study->{columns};
  # the sort is important for the signature
  my $this_entity_columns = [ sort grep { $column_config->{$_}{describes} eq $entity->{name} } @$column_keys ];
  my $signature = join '::', map { $row->{$_} // '' } @$this_entity_columns;

  return '' if ($signature =~ /^:*$/); # if all the columns were empty the signature will be empty or just colons

  # Prefix the signature by the parent entity ID.
  # We don't need to include ALL parent entity IDs back to the root because
  # the direct parent's ID should already be unique.
  my $parent_entity = $entity_to_parent->{$entity};
  if (defined $parent_entity) { # root entity has no parent
    my $parent_id_column = find_id_column_name($column_config, $parent_entity);
    my $parent_id = $row->{$parent_id_column};
    unless (defined $parent_id) {
      my $row_snippet = encode_json({ map { $_ => $row->{$_} } @$this_entity_columns });
      die "Can't find a parent ID for make_auto_entity_id in column >$parent_id_column< at row: $row_snippet\n";
    }
    $signature = join '::', $parent_id, $signature;
  }

  if (!$seen_signatures{$entity->{name}}{$signature}) {
    my $new_id = sprintf '%s-%05d', $entity->{name}, 1 + keys %{$seen_signatures{$entity->{name}}};
    $new_id =~ s/ /_/g; # change spaces to underscores
    $seen_signatures{$entity->{name}}{$signature} = $new_id;
  }
  return $seen_signatures{$entity->{name}}{$signature};
}


# just return the prefix before the underscore or 'TERM' as fallback
sub source_ref {
  my ($term_id) = @_;
  my ($prefix, $acc) = split /_/, $term_id;
  return $prefix || 'TERM';
}


sub find_id_column_name {
  my ($column_config, $entity) = @_;

  # Find the column name that matches the criteria
  my ($id_column_name) = grep {
    !$column_config->{$_}{ignore} &&
      $column_config->{$_}{describes} eq $entity->{name} &&
      $column_config->{$_}{value_type} eq 'id'
    } keys %$column_config;

  return $id_column_name;
}

sub ordered_hashref {
  my $ref = {};
  tie %{$ref}, 'Tie::Hash::Indexed';
  return $ref;
}

# make a hashref lookup for all term IDs matching "http://purl.obolibrary.org/obo/ABC_001234"
sub get_valid_term_ids {
  my @urls = @_;
  my $lookup = {};
  warn "Scanning these OWL files for valid term IDs:\n".join("\n", @urls)."\n";
  foreach my $url (@urls) {
    my $file_contents = get($url);
    while ($file_contents =~ m{http://purl.obolibrary.org/obo/(\w+_\d+)}g) {
      $lookup->{$1}++;
    }
  }
  my $n = scalar %$lookup;
  warn "Done. Found $n IDs.\n";
  return $lookup;
}


sub grep_keys_get_values {
  my ($data, $pattern) = @_;
  my @result;
  if (ref $data eq 'HASH') {
    for my $key (keys %$data) {
      my $value = $data->{$key};
      push @result, $value if ($key =~ /$pattern/);
      push @result, grep_keys_get_values($value, $pattern);
    }
  } elsif (ref $data eq 'ARRAY') {
    for my $element (@$data) {
      push @result, grep_keys_get_values($element, $pattern);
    }
  }
  return @result;
}
