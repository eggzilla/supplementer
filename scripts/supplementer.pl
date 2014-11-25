#!/usr/bin/perl

use version; our $VERSION = qv('0.01');
use strict;
use warnings;
use Template;
use Cwd;
use File::Basename;
use IPC::Cmd qw[can_run run run_forked];

sub make_supplements{
  my ($csv_file_path, $goi_root, $html_destination_path, $base_URL, $log_path) = @_;
  #check arguments
  die ("ERROR $fasta_file_path does not exist\n") unless (-e $fasta_file_path);
  die ("ERROR $html_destination_path does not exist\n") unless (-d $assembly_hub_destination_path);
  die ("ERROR no URL (network location) provided" unless(defined $base_URL);
  die ("ERROR $log_path does not exist\n") unless (-e $log_path);
  #ensure that base_URL ends with slash
  $base_URL =~ s!/*$!/!;  

  #check program dependencies


  #create html directory structure

  #template definition
  my $template = Template->new({
                                INCLUDE_PATH => ["$template_path"],
                                RELATIVE=>1,
  });

  #construct index.hmtl
  my $index_path = $html_destination_path. "/index.html";
  my $index_file = 'index.html';
  my $index_vars =
    {
     foo => "$foo",
    };
  $template->process($index_file,$index_vars,$index_path) || die "Template process failed: ", $template->error(), "\n";

  #construct gene of interest goi.html
  my $goi_path = $assembly_hub_directory . "/goi.html";
  my $goi_file = 'goi.html';
  my $goi_vars =
    {
     bar => "$bar"
    };
  $template->process($goi_file,$goi_vars,$goi_path) || die "Template process failed: ", $template->error(), "\n";
}
