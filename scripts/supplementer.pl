#!/usr/bin/perl

#Script supplementer.pl;
#Last changed Time-stamp: <2014-11-26 11:53:50 fall> by joerg
###############
###Use stuff
###############
use version; our $VERSION = qv('0.01');
use strict;
use warnings;
use Template;
use Cwd;
use File::Basename;
use File::Path qw(make_path remove_tree);
use IPC::Cmd qw[can_run run run_forked];
use PerlIO::gzip;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Cwd;
use Text::Iconv;
 # my $converter = Text::Iconv -> new ("utf-8", "windows-1251");
 # Text::Iconv is not really required.
 # This can be any object with the convert method. Or nothing.
use Spreadsheet::XLSX;
use DateTime qw();

###############
###Variables
###############
my $VERBOSE = 0;
my ( $dir, $odir, $files, $outfiles );

###############
###Command Line Options
###############

pod2usage(-verbose => 0)
    unless GetOptions(
	"dir|d=s"      => \$dir,
	"odir|o=s"     => \$odir,
	"files|f=s"    => \$files,
	"ofiles|t=s"   => \$outfiles,
	"help|h"       => sub{pod2usage(-verbose => 1)},
	"man|m"	       => sub{pod2usage(-verbose => 2)},
	"verbose"      => sub{ $VERBOSE++ }
    );

$dir  =	 cwd() unless ($dir);
my @csv = split(/\,/,$files);
my $today = DateTime->now->strftime('%d%m%Y');
$odir =	 "$dir"."\/Supplements_$today/" unless $odir;
$dir  =~ s/ //g;
$odir =~ s/ //g;
($dir) or die "No working directory chosen!\n";

my $pid = $$;
(my $job = `cat /proc/$pid/cmdline`)=~ s/\0/ /g;
print STDERR $job,"\n";

###############
###Main Stuff
###############

chdir ($dir) or die "Directory $dir could not be found!\n";

if (!-d $odir){
    print STDERR "Creating output directory $odir!\n";
    make_path ($odir) or die "Error creating directory: $odir";
}

print STDERR "Reading input!\n";

foreach my $file (@csv){
    die "File $file could not be openend!\n" unless (-e $file);
    #### Read fields for html from file

    my ($goi, @synonyms, @pathways); 

    if ($file =~ /.csv/){
	print STDERR "Parsing csv!\n";
	open (LIST,"<:gzip(autopop)","$file");
	while(<LIST>){
	    chomp(my $line = $_);
	    my @fields  = split(/,/,$line);
	    $goi	= $fields[0];
	    @synonyms	= split(",",$fields[1]);
	    @pathways	= split(",",$fields[2]);
	    
 
	}
    }
    elsif ($file =~ /.xlsx/){
	print STDERR "Parsing Excel sheet!\n";
	my $excel = Spreadsheet::XLSX -> new ($file);
	foreach my $sheet (@{$excel -> {Worksheet}}) {	    
	    next unless ($sheet->{Name} eq 'APG' || $sheet->{Name} eq 'GOI');
	    printf("Sheet: %s\n", $sheet->{Name});    
	    $sheet -> {MaxRow} ||= $sheet -> {MinRow};	    
	    foreach my $row ( $sheet -> {MinRow} .. $sheet -> {MaxRow}) {
		print STDERR "ROW: $row\n";
               $sheet -> {MaxCol} ||= $sheet -> {MinCol};                
                foreach my $col ($sheet -> {MinCol} ..  $sheet -> {MaxCol}) {		    
		    print STDERR "COL: $col\n";
		    my $cell = $sheet -> {Cells} [$row] [$col];		    
		    if ($cell) {
			printf("( %s , %s ) => %s\n", $row, $col, $cell -> {Val});
		    }	    
                }		
	    }	    
	}
    }
}

#sub make_supplements{
#  my ($csv_file_path, $goi_root, $html_destination_path, $base_URL, $log_path) = @_;
#  #check arguments
#  die ("ERROR $fasta_file_path does not exist\n") unless (-e $fasta_file_path);
#  die ("ERROR $html_destination_path does not exist\n") unless (-d $assembly_hub_destination_path);
#  die ("ERROR no URL (network location) provided" unless(defined $base_URL);
#  die ("ERROR $log_path does not exist\n") unless (-e $log_path);
#  #ensure that base_URL ends with slash
#  $base_URL =~ s!/*$!/!;  
#
#  #check program dependencies
#
#
#  #create html directory structure
#
#  #template definition
#  my $template = Template->new({
#                                INCLUDE_PATH => ["$template_path"],
#                                RELATIVE=>1,
#  });
#
#  #construct index.hmtl
#  my $index_path = $html_destination_path. "/index.html";
#  my $index_file = 'index.html';
#  my $index_vars =
#    {
#     foo => "$foo",
#    };
#  $template->process($index_file,$index_vars,$index_path) || die "Template process failed: ", $template->error(), "\n";
#
#  #construct gene of interest goi.html
#  my $goi_path = $assembly_hub_directory . "/goi.html";
#  my $goi_file = 'goi.html';
#  my $goi_vars =
#    {
#     bar => "$bar"
#    };
#  $template->process($goi_file,$goi_vars,$goi_path) || die "Template process failed: ", $template->error(), "\n";
#}
