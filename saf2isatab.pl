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

GetOptions("output_directory|output-directory=s" => \$output_dir

	  );

my ($config_filename, $saf_filename) = @ARGV;

die "must provide two filenames on commandline: config_file saf_data_file\n"
    unless ($config_filename && $saf_filename && -e $config_filename && -e $saf_filename);

my $userConfig = LoadFile($config_filename);
die "problem reading '$config_filename'\n" unless $userConfig;

# merge the user config into the default config

my $config = merge $defaultConfig, $userConfig;


my $hashified = Text::CSV::Hashify->new({
					 file => $saf_filename,
					 key => 'sample_ID',
					 sep_char => $saf_filename =~ /\.csv$/ ? ',' : "\t",
					});

my $sample_IDs = $hashified->keys;
my $columns = $hashified->fields;

print "@$columns\n";

my $isa_writer = Bio::Parser::ISATab->new(directory => $output_dir, protocols_first=>1);
# need to set $config->{study_file_name} somewhere to the s_samples file
$isa_writer->write( { ontologies => [], studies => [ $config ] } );

