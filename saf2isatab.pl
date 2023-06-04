#!/usr/bin/env perl
# -*- mode: cperl -*- 


#
# usage: ./saf2isatab.pl config.yaml saf-data.tsv
#
# (can handle tsv or csv data)
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
my $column_config = delete $config->{columns};
# print Dumper($column_config); exit;
# and some other non ISA-Tab bits
my $study_species = delete $config->{study_species};
my $study_developmental_stages = delete $config->{study_developmental_stages};
my $study_sexes = delete $config->{study_sexes};
my $study_terms = delete $config->{study_terms};
my $location_qualifiers = delete $config->{location_qualifiers};

# what remains is ISA-Tab Study material
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
validate_column_config($column_config, $flat_entities);


# load the actual data
my $hashified = Text::CSV::Hashify->new({
					 file => $saf_filename,
					 key => 'sample_ID',
					 sep_char => $saf_filename =~ /\.csv$/ ? ',' : "\t",
					});

my $sample_IDs = $hashified->keys;
my $column_keys = $hashified->fields;

# make sure no required columns are missing from column_keys
# and warn about any unconfigured columns
validate_columns($column_keys, $column_config);

$study->{study_file_name} = 's_samples.txt';

my $sources = $study->{sources} = {};
my $study_assays = $study->{study_assays} = [];
my $study_protocols = $study->{study_protocols} //= [];

foreach my $sample_ID (@$sample_IDs) {
  warn "trying >$sample_ID<\n";
  my $row = $hashified->record($sample_ID);

  # add material entities (descending the entity graph into assay entities also)
  add_material($root_entity, $row, $sources, $study_assays, $study_protocols);
}


#
# write the ISA-Tab!
#
my $isa_writer = Bio::Parser::ISATab->new(directory => $output_dir, protocols_first=>1);
$isa_writer->write( { ontologies => [], studies => [ $study ] } );



sub add_material {
  my ($entity, $row, $isaref, $study_assays, $study_protocols) = @_;

  # figure out an ID for this entity
  # this is simple if all *_ID fields are mandatory
  # otherwise we'll need to generate default IDs (which we can add later if providing, e.g. location_ID is a hassle)

  my ($id_column_name) =
    grep { $column_config->{$_}{describes} eq $entity->{name} &&
	   $column_config->{$_}{value_type} eq 'id' }
    keys %$column_config;

  die "FATAL ERROR: couldn't find ID column name for $entity->{name}\n"
    unless ($id_column_name);
  
  warn "id column is $id_column_name\n";


  # recurse down entity tree
  foreach my $child_entity (@{$entity->{children}}) {
    # check material or assay
    

  }

}


sub validate_column_config {
  my ($column_config, $flat_entities) = @_;

  # check that every column definition has `describes` and `value_type`

  my @bad = grep { !$column_config->{$_}{describes} || !$column_config->{$_}{value_type} } keys %$column_config;

  die "FATAL ERROR: column configuration is missing required attributes 'describes' and/or 'value_type' for columns: ".join(', ', @bad)."\n" if (@bad);

  # now check that all the columns `describe` an existing entity
  my %entity_names = map { ($_->{name} => 1) } @$flat_entities;
  my @worse = grep { !$entity_names{$column_config->{$_}{describes}} } keys %$column_config;

  die "FATAL ERROR: these columns 'describe' entities that do not exist: ".join(', ', @worse)."\n"
    if (@worse);

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


sub flatten {
  my ($entity) = @_;
  return [ $entity, map { @{ flatten($_) } } @{$entity->{children}} ];
}
