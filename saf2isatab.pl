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
# and some other non ISA-Tab bits
my $study_species = delete $config->{study_species};
my $study_developmental_stages = delete $config->{study_developmental_stages};
my $study_sexes = delete $config->{study_sexes};
my $study_terms = delete $config->{study_terms};

# what remains is ISA-Tab Study material
my $study = $config;

# load the entity graph
my $entities = LoadFile($entities_filename);

die "root entity (in '$entities_filename') must be a material type entity\n"
  unless ($entities->[0]{type} eq 'material');

# load the actual data
my $hashified = Text::CSV::Hashify->new({
					 file => $saf_filename,
					 key => 'sample_ID',
					 sep_char => $saf_filename =~ /\.csv$/ ? ',' : "\t",
					});

my $sample_IDs = $hashified->keys;
my $column_keys = $hashified->fields;

print "@$column_keys\n";

$study->{study_file_name} = 's_samples.txt';

my $sources = $study->{sources} = {};
my $study_assays = $study->{study_assays} = [];
my $study_protocols = $study->{study_protocols} //= [];

foreach my $sample_ID (@$sample_IDs) {
  warn "trying >$sample_ID<\n";
  my $row = $hashified->record($sample_ID);

}


#
# write the ISA-Tab!
#
my $isa_writer = Bio::Parser::ISATab->new(directory => $output_dir, protocols_first=>1);
$isa_writer->write( { ontologies => [], studies => [ $study ] } );

