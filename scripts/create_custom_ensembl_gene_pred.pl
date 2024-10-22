#!/usr/bin/env perl

#To use this script you will need to figure out which GENCODE version corresponds to a particular version of Ensembl
#e.g. 
#   For Ensembl v105, visit the 105 archive site: https://dec2021.archive.ensembl.org/Homo_sapiens/Info/Annotation
#   Which tells you that "Database version 105.38" corresponds to "Gencode version GENCODE 39"

#To obtain the starting material for a genePred file, use the UCSC Genome Table Browser.
# https://genome.ucsc.edu/cgi-bin/hgTables
# Make sure the following options are selected (using GENCODE 39 as an example) and then hit "Get output"
#   Genome:Human 
#   Assembly:GRCh38
#   Group:Genes and Gene Predictions 
#   Track:All GENCODE V39
#   Table: Basic
#   Region:Genome
#   Output format:Selected fields from primary and related tables
#   Output field separator: tab
#   File type returned: Gzip compressed

#On the next page, link the tables: wgEncodeGencodeAttrsV39, wgEncodeGencodeTagV39, wgEncodeGencodeTranscriptionSupportLevelV39
#And select the following fields:
#   wgEncodeGencodeAttrsV39.geneId, wgEncodeGencodeAttrsV39.transcriptType, 
#   wgEncodeGencodeTagV39.tag
#   wgEncodeGencodeTranscriptionSupportLevelV39.level

#Note you may need Mozilla::CA to get access to the https request to work
#Or you might try this instead: LWP::Protocol::https

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

my $ucsc_pred_file = '';
my $custom_pred_file = '';
my $build = '';
my $gencode_version = '';
my $protein_coding = '';
my $tsl1 = '';

GetOptions ('ucsc_pred_file=s'=>\$ucsc_pred_file, 'custom_pred_file=s'=>\$custom_pred_file,
            'build=s'=>\$build, 'gencode_version=s'=>\$gencode_version, 'protein_coding=i'=>\$protein_coding,
            'tsl1=i'=>\$tsl1);

my $usage=<<INFO;
  Example usage: 
  
  ./createCustomEnsemblGenePred.pl --ucsc_pred_file='Ensembl105_GRCh38_UcscGenePred.ensGene.gz' --custom_pred_file='Ensembl105_GRCh38_UcscGenePred_Custom_Coding.ensGene' --build='hg38' --gencode_version='V39' --protein_coding=1

  --ucsc_pred_file              Gene annotation file for ensembl transcripts obtained from UCSC table browser
                                Supplemented with additional columns to allow filtering to neoantigen relevant transcript models
                                See comments in this script for details

  --custom_pred_file            Gene prediction file, compatible with IGV for use in genomic review of neoantigen candidates

  --build                       Genome build name used in UCSC table headers (default = hg38)

  --gencode_version             GENCODE version name used in UCSC table headers (default = V39)

  --protein_coding              Limit to protein coding transcripts (0 or 1, default = 0)

  --tsl1                        Limit to TSL=1 transcripts (0 or 1, default = 0)

INFO

unless (-e $ucsc_pred_file && $custom_pred_file){
  print RED, "\n\nRequired parameter missing or invalid", RESET;
  print GREEN, "\n\n$usage", RESET;
  exit(1);
}
unless ($build){ $build = "hg38";}
unless ($gencode_version){ $gencode_version = "V39";}
unless ($protein_coding){ $protein_coding = 0;}
unless ($tsl1){ $tsl1 = 0;}


#Create a new gene pred file based on a UCSC table browser query to include only the genes meeting user specified criteria
#Modify the gene pred file to use names that more useful for review: (Ensembl transcript IDs with other info appended)
print GREEN, "\n\nCreating a customized gene pred file for Ensembl transcripts", RESET;

if($protein_coding){
  print YELLOW, "\nLimiting to protein coding transcripts", RESET;
}

open(OUTPRED, ">$custom_pred_file") || die RED, "\n\nCould not open output gene pred file: $custom_pred_file\n\n";
open(INPRED, "gunzip -c $ucsc_pred_file |") || die RED, "\n\nCould not open input gene pred file: $ucsc_pred_file\n\n", RESET;

#Define expected column names
my %columns;
my $enst_col_name = "$build".".wgEncodeGencodeBasic"."$gencode_version".".name";
my $gene_symbol_col_name = "$build".".wgEncodeGencodeBasic"."$gencode_version".".name2";
my $ensg_id_col_name = "$build".".wgEncodeGencodeAttrs"."$gencode_version".".geneId";
my $transcript_type_col_name = "$build".".wgEncodeGencodeAttrs"."$gencode_version".".transcriptType";
my $transcript_tags_col_name = "$build".".wgEncodeGencodeTag"."$gencode_version".".tag";
my $tsl_col_name = "$build".".wgEncodeGencodeTranscriptionSupportLevel"."$gencode_version".".level";
my @expected_columns = ($enst_col_name, $gene_symbol_col_name, $ensg_id_col_name, $transcript_type_col_name, $transcript_tags_col_name, $tsl_col_name);

while(<INPRED>){
  chomp $_;
  my @line = split("\t", $_);
  my @pred_line = @line[0..15];

  if ($_ =~ /^\#/){
    my $string = join("\t", @pred_line);
    print OUTPRED "$string\n";

    my $pos = 0;
    foreach my $head (@line){
      $columns{$head}{pos} = $pos;
      $pos++;
    }
    foreach my $expected_col_name (@expected_columns){
      unless (defined $columns{$expected_col_name}){
        print RED, "\n\nCould not find expected column in header of input file: $expected_col_name. Check input file and build/gencode_version options\n\n", RESET; 
        exit;
      }
    }
    next;
  }
  
  my $enst = $line[$columns{$enst_col_name}{pos}];
  my $enst_compact = $enst;
  $enst_compact =~ s/ENST0+/ENST~/;
  my $symbol = $line[$columns{$gene_symbol_col_name}{pos}];
  my $ensg_id = $line[$columns{$ensg_id_col_name}{pos}];
  my $transcript_type = $line[$columns{$transcript_type_col_name}{pos}];
  my $transcript_tags = $line[$columns{$transcript_tags_col_name}{pos}];
  my $tsl = $line[$columns{$tsl_col_name}{pos}];
  
  #If requested by the user, limit to "protein_coding" transcripts
  if ($protein_coding){
    next unless ($transcript_type eq "protein_coding");
  }

  #If requested by the user, limit to TSL level=1 transcripts
  if ($tsl1){
    next unless ($tsl eq "1");
  }

  #Append a marking to the Ensembl MANE_select transcript where possible
  my $marker="";
  if ($transcript_tags =~ /MANE_Select/i){
     $marker = "***";
  }

  #If requested by the user, limit to transcripts corresponding to the BEST transcript from a pVACseq aggregate report file


  #Append ensembl version info to the enst id
  #Also create a more useful gene name that shows: symbol, ENST ID, TSL
  my $custom_name = $marker . $symbol . "_" . $enst_compact . " (TSL=$tsl)" . $marker;
  $pred_line[12] = $custom_name;

  my $string = join("\t", @pred_line);
  print OUTPRED "$string\n";
}
close(INPRED);
close(OUTPRED);

print GREEN, "\n\nNew gene pred file created: $custom_pred_file\n\n", RESET;

exit;

