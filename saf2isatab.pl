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
#- do range checking on lat/long and date format checks
#
#- collect and report all errors rather than die for each one
#
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
my $default_config_filename = $FindBin::Bin."/default-column-config.yaml";
my $defaultConfig = LoadFile($default_config_filename);

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
validate_columns($column_keys, $column_config);

# append mandatory columns that have default values
$column_keys = add_mandatory_columns($column_keys, $column_config);


$study->{study_file_name} = 's_samples.txt';

my $study_assays = $study->{study_assays} = [];
my $study_protocols = $study->{study_protocols} //= [];

foreach my $sample_ID (@$sample_IDs) {
  warn "trying >$sample_ID<\n";
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
    grep { $column_config->{$_}{describes} eq $entity->{name} &&
	   $column_config->{$_}{value_type} eq 'id' }
    keys %$column_config;

  die "FATAL ERROR: couldn't find ID column name for $entity->{name}\n"
    unless ($id_column_name);

  my $entity_id = $row->{$id_column_name};
  # warn "id column is $id_column_name and got $entity_id\n";


  # make a hashref to put the new material nodes in, e.g.
  # $study->{sources}{some_source_id} = ...

  $isaref->{$entity->{isa_key}} //= {};
  my $new_isaref = $isaref->{$entity->{isa_key}}{$entity_id};

  # only make a node once (i.e. don't reprocess a location multiple times for each collection and sample)
  if (!defined $new_isaref) {
    $new_isaref = $isaref->{$entity->{isa_key}}{$entity_id} //= {};

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
	my $characteristics = $new_isaref->{characteristics}{"$column (TERM:$col_term)"} //= {};

	# if it's a plain text/number/date value then it's a simple case
	if ($col_config->{value_type} =~ /^(number|string|date|latitude|longitude)$/) {
	  $characteristics->{value} = $value;

	  # ontology term values require a lookup from text to term:
	} elsif ($col_config->{value_type} eq 'term') {
	  # get the lookup hash (already validated - no need to check success)
	  my $lookup = $config->{$col_config->{term_lookup} // 'study_terms'};
	  my $value_term_id = $lookup->{$value};
	  if ($value_term_id) {
	    $characteristics->{value} = $value;
	    $characteristics->{term_source_ref} = 'TERM';
	    $characteristics->{term_accession_number} = $value_term_id;
	  } else {
	    die sprintf "FATAL ERROR: value '%s' not found in '%s' term lookup\n", $value, $col_config->{term_lookup} // 'study_terms';
	  }
	}
      }
    }
  }

  # recurse down entity tree
  foreach my $child_entity (@{$entity->{children}}) {
    # check material or assay
    if ($child_entity->{type} eq 'material') {
      add_material($child_entity, $row, $new_isaref, $column_keys, $config_and_study);
    }
  }

}


sub validate_config {
  my ($config, $flat_entities) = @_;

  my $column_config = $config->{columns};
  # check that every column definition has `describes` and `value_type`
  my @bad = grep { !$column_config->{$_}{describes} || !$column_config->{$_}{value_type} } keys %$column_config;
  die "FATAL ERROR: column configuration is missing required attributes 'describes' and/or 'value_type' for columns: ".join(', ', @bad)."\n" if (@bad);

  # now check that all the columns `describe` an existing entity
  my %entity_names = map { ($_->{name} => 1) } @$flat_entities;
  my @worse = grep { !$entity_names{$column_config->{$_}{describes}} } keys %$column_config;
  die "FATAL ERROR: these columns 'describe' entities that do not exist: ".join(', ', @worse)."\n"
    if (@worse);

  # check that every non-deprecated term|number|string|date|latitude|longitude column has a column_term
  my @terrible = grep {
    $column_config->{$_}{value_type} =~ /^(term|number|string|date|latitude|longitude)$/ &&
    !$column_config->{$_}{deprecated} &&
    !$column_config->{$_}{column_term}
  } keys %$column_config;
  die "FATAL ERROR: these columns don't have a column_term (e.g. a variable IRI)".join(', ', @terrible)."\n"
    if (@terrible);

  # check that any term_lookup values exist in the $config hash as first level keys
  my @dreadful = grep {
    $column_config->{$_}{term_lookup} && !exists $config->{$column_config->{$_}{term_lookup}}
  } keys %$column_config;
  die "FATAL ERROR: these columns have term_lookup values that are not defined in the config file:".join(', ', @dreadful)."\n"
    if (@dreadful);

  # this has side effects!
  # add `required: true` to any column that doesn't have it
  map { $_->{required} //= 1 } values %$column_config;
}


sub validate_columns {
  my ($column_keys, $column_config) = @_;

  # make a hash look-up
  my %column_keys = map { ($_ => 1) } @$column_keys;

  # are all the required columns in the input file's header?
  #
  my @missing = grep {
      ($column_config->{$_}{required} && !defined $column_config->{$_}{default}) && !$column_keys{$_}
  } keys %$column_config;

  die "FATAL ERROR: the following required columns are missing from the input file: ".join(', ', @missing)."\n"
    if (@missing);


  my @unconfigured = grep {
    !exists $column_config->{$_}
  } @$column_keys;

  warn "WARNING: the following data file columns are not configured and will be ignored: ".join(', ', @unconfigured)."\n"
    if (@unconfigured);
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
      $column_config->{$_}{required} && # and it's required
      defined $column_config->{$_}{default} # and there's a default provided
    } keys %$column_config;

  return $result;
}

sub flatten {
  my ($entity) = @_;
  return [ $entity, map { @{ flatten($_) } } @{$entity->{children}} ];
}
