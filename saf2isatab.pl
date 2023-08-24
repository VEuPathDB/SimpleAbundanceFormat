#!/usr/bin/env perl
# -*- mode: cperl -*- 


#
# usage: ./saf2isatab.pl config.yaml saf-data.tsv
#
# (can handle tsv or csv data)
#

#
# TO DO:
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
#
# LIMITATIONS:
#
#- multiple protocols per assay/material processing aren't handled, but hopefully we won't need them
#



use strict;
use warnings;

use Getopt::Long;
use YAML::XS qw/LoadFile/;
use Text::CSV::Hashify; # previously we appended a version number 0.11
use Bio::Parser::ISATab;
use FindBin;
use Hash::Merge::Simple qw/merge/;
use Data::Dumper;

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

GetOptions("output_directory|output-directory=s" => \$output_dir,
	   "entities=s" => \$entities_filename,
	  );

my ($config_filename, $saf_filename) = @ARGV;

die "must provide two filenames on commandline: config_file saf_data_file\n"
    unless ($config_filename && $saf_filename && -e $config_filename && -e $saf_filename);

my $userConfig = LoadFile($config_filename);
die "problem reading '$config_filename'\n" unless $userConfig;

# merge the user config into the default config
my $config = merge $defaultConfig, $userConfig;

# pull out the column config
my $column_config = $config->{columns};
# print Dumper($column_config); exit;

# a subset of the config is also the basis of the ISA-Tab study datastructure
# the ISA-Tab writer ignores the extra config information
my $study = $config;

# the YAML merge can't handle arrays, so we append the default_study_protocols onto the main protocols
# (with no checks for duplicates...)
push @{$study->{study_protocols}}, @{$config->{default_study_protocols}};

# load the entity graph
my $entities = LoadFile($entities_filename);

die "FATAL ERROR: there must be exactly one root entity (in '$entities_filename')\n"
  if (@$entities != 1);
die "FATAL ERROR: root entity (in '$entities_filename') must be a material type entity\n"
  unless ($entities->[0]{type} eq 'material');

my $root_entity = $entities->[0];
my $flat_entities = flatten($root_entity);

# make sure column config contains required info
validate_config($config, $flat_entities);


# load the actual data
my $hashified = Text::CSV::Hashify->new({
					 file => $saf_filename,
					 key => 'sample_ID',
					 sep_char => $saf_filename =~ /\.csv$/ ? ',' : "\t",
					});

# an arrayref of the sample IDs in the input file
my $sample_IDs = $hashified->keys;

# an arrayref of the column headings in the input file
my $column_keys = $hashified->fields;

# make sure no required columns are missing from column_keys
# and warn about any unconfigured columns
$column_keys = validate_columns($column_keys, $column_config);

# append mandatory columns that have default values
$column_keys = add_mandatory_columns($column_keys, $column_config);


$study->{study_file_name} = 's_samples.txt';

my $study_assays = $study->{study_assays} = [];
my $study_protocols = $study->{study_protocols} //= [];

foreach my $sample_ID (@$sample_IDs) {
  warn "processing sample: $sample_ID\n";
  my $row = $hashified->record($sample_ID);

  # add material entities (descending the entity graph into assay entities also)
  add_material($root_entity, $row, $study, $column_keys, $config);
}


#
# write the ISA-Tab!
#
my $isa_writer = Bio::Parser::ISATab->new(directory => $output_dir, protocols_first=>1);
$isa_writer->write( { ontologies => [], studies => [ $study ] } );



sub add_material {
  my ($entity, $row, $isaref, $column_keys, $config_and_study) = @_;

  my $column_config = $config_and_study->{columns};
  # figure out an ID for this entity
  # this is simple if all *_ID fields are mandatory
  # otherwise we'll need to generate default IDs (which we can add later if providing, e.g. location_ID is a hassle)

  my ($id_column_name) =
    grep { !$column_config->{$_}{ignore} &&
	   $column_config->{$_}{describes} eq $entity->{name} &&
	   $column_config->{$_}{value_type} eq 'id' }
    keys %$column_config;

  die "FATAL ERROR: couldn't find ID column name for $entity->{name}\n"
    unless ($id_column_name);

  my $entity_id = $row->{$id_column_name};
  # warn "id column is $id_column_name and got $entity_id\n";

  if (!defined $entity_id && $column_config->{$id_column_name}{default} eq '__AUTO__') {
    $entity_id = make_auto_entity_id($entity, $row, $column_keys, $config_and_study);
  }

  # make a hashref to put the new material nodes in, e.g.
  # $study->{sources}{some_source_id} = ...
  $isaref->{$entity->{isa_key}} //= {}; # could/should be ordered hashref
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

    # handle the characteristics/variables (not comments or protocols)
    if ($col_term) {
      my $characteristics = $isaref->{characteristics}{"$column ($col_term)"} //= {};

      # if it's a plain text/number/date value then it's a simple case
      # multivalued values can be left as they are
      if ($col_config->{value_type} =~ /^(number|string|date|latitude|longitude)$/) {
	$characteristics->{value} = $value;

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
	  } else {
	    die sprintf "FATAL ERROR: value '%s' not found in '%s' term lookup\n", $value, $col_config->{term_lookup} // 'study_terms';
	  }
	}
	$characteristics->{value} = join $default_isatab_delimiter, @values;
	$characteristics->{term_source_ref} = join $default_isatab_delimiter, @term_source_refs;
	$characteristics->{term_accession_number} = join $default_isatab_delimiter, @term_accession_numbers;
      }
    } elsif ($col_config->{value_type} eq 'protocol_ref') {
      # check that the protocol ref in $value is in the study_protocols
      my @ok = grep { $_->{study_protocol_name} eq $value } @{$config_and_study->{study_protocols}};
      if (@ok) {
	$isaref->{protocols}{$value} = {};
      } else {
	die "FATAL ERROR: protocol ref '$value' not found in study_protocols\n";
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

    # now add the link from sample ID to the actual assay
    my $assay_id = $row->{sample_ID}.'.'.($row_protocol_ref || $column_protocol_ref);
    my $assay = $study_assay->{samples}{$row->{sample_ID}}{assays}{$assay_id} //= {};
    # add the column-wise protocol ref explicitly
    $assay->{protocols}{$column_protocol_ref} = {}
      if ($column_protocol_ref);
    # otherwise row-wise protocol ref will be added from the input $row by add_column_data
    add_column_data($assay, $column_keys, $row, $entity, $config_and_study);
  }
}

sub validate_config {
  my ($config, $flat_entities) = @_;

  my $column_config = $config->{columns};
  # check that every column definition has `describes` and `value_type`
  my @bad = grep { (!$column_config->{$_}{describes} || !$column_config->{$_}{value_type}) && !$column_config->{$_}{ignore} } keys %$column_config;
  die "FATAL ERROR: column configuration is missing required attributes 'describes' and/or 'value_type' for columns: ".join(', ', @bad)."\n" if (@bad);

  # now check that all the columns `describe` an existing entity
  my %entity_names = map { ($_->{name} => 1) } @$flat_entities;
  my @worse = grep { !$column_config->{$_}{ignore} && !$entity_names{$column_config->{$_}{describes}} } keys %$column_config;
  die "FATAL ERROR: these columns 'describe' entities that do not exist: ".join(', ', @worse)."\n"
    if (@worse);

  # check that every non-deprecated term|number|string|date|latitude|longitude column has a column_term
  my @terrible = grep {
    !$column_config->{$_}{ignore} &&
    $column_config->{$_}{value_type} =~ /^(term|number|string|date|latitude|longitude)$/ &&
    !$column_config->{$_}{deprecated} &&
    !$column_config->{$_}{column_term}
  } keys %$column_config;
  die "FATAL ERROR: these columns don't have a column_term (e.g. a variable IRI): ".join(', ', @terrible)."\n"
    if (@terrible);

  # check that any term_lookup values exist in the $config hash as first level keys
  my @dreadful = grep {
    !$column_config->{$_}{ignore} &&
    $column_config->{$_}{term_lookup} &&
    !exists $config->{$column_config->{$_}{term_lookup}}
  } keys %$column_config;
  die "FATAL ERROR: these columns have term_lookup values that are not defined in the config file:".join(', ', @dreadful)."\n"
    if (@dreadful);


  # check that any `protocol` values for assay variables are in the study_protocols
  # (nested grep is not too pretty)
  my @awful = grep {
    my $column = $_;
    $column_config->{$column}{protocol} &&
    !grep { $_->{study_protocol_name} eq $column_config->{$column}{protocol} } @{$config->{study_protocols}}
  } keys %$column_config;
  die "FATAL ERROR: these columns contain `protocol` values that are not in study_protocols: ".join(', ', @awful)."\n"
    if (@awful);

  ### the following have side effects!

  # add the default `required: true` to any column that doesn't have it
  map { $_->{required} //= 1 } values %$column_config;

  # add default delimiter for multivalued variables
  map { $_->{delimiter} //= $default_input_delimiter } grep { $_->{multivalued} } values %$column_config;
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
     samples => {},
    };

  push @{$config_and_study->{study_assays}}, $study_assay;
  return $study_assay;
}


#
# generate an entity_id based on the unique signature of all its column values
#

# we use the size of this hash to generate the entity ID serial number
my %seen_signatures;

sub make_auto_entity_id {
  my ($entity, $row, $column_keys, $config_and_study) = @_;

  my $column_config = $config_and_study->{columns};
  # the sort is important for the signature
  my $relevant_columns = [ sort grep { $column_config->{$_}{describes} eq $entity->{name} } @$column_keys ];
  my $signature = join '::', map { $row->{$_} // '' } @$relevant_columns;

  if (!$seen_signatures{$signature}) {
    my $new_id = sprintf '%s-%05d', $entity->{name}, 1 + keys %seen_signatures;
    $seen_signatures{$signature} = $new_id;
  }
  return $seen_signatures{$signature};
}


# just return the prefix before the underscore or 'TERM' as fallback
sub source_ref {
  my ($term_id) = @_;
  my ($prefix, $acc) = split /_/, $term_id;
  return $prefix || 'TERM';
}
