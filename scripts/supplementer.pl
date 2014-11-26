#!/usr/bin/perl

#Script supplementer.pl;
#Last changed Time-stamp: <2014-11-26 13:02:29 fall> by joerg
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
use Spreadsheet::ParseXLSX;
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

    my ($wdir,$filetoparse) = split("\/",$file,2);
#    print STDERR "$wdir, $filetoparse\n";
    chdir ($wdir) or die $!;
    my %entries; 

    if ($filetoparse =~ /.csv/){
	print STDERR "Parsing csv $filetoparse!\n";
	open (LIST,"<:gzip(autopop)","$file");
	while(<LIST>){
	    next unless ($_ =~ /.goi./ || $_ =~ /.apk./);
	    chomp(my $line	  = $_);
	    my @fields		  = split(/\t/,$line);
	    my $goi		  = $fields[0];
	    my $hacker		  = $fields[1];
	    my $gene		  = $fields[3];
	    my $duplicate	  = $fields[4];   
	    my @synonyms	  = split(",",$fields[5]);
	    my @pathways	  = split(",",$fields[6]);
	    my @literature	  = split(",",$fields[7]);
	    my $igvs		  = $fields[8];
	    my $mock		  = $fields[9];
 	    my $ebola		  = $fields[10];
	    my $marburg		  = $fields[11];
	    my $profile		  = $fields[12];
	    my $homolog		  = $fields[13];
	    my $raemock		  = $fields[14];
 	    my $raeebola	  = $fields[15];
	    my $raemarburg	  = $fields[16];
	    my $raeprofile	  = $fields[17];
	    my $ucscs		  = $fields[18];
	    my $ucsc_conservation = $fields[19];
	    my $hg_SNPs		  = $fields[20];
	    my $rae_SNPs	  = $fields[21];
	    my $sashimi		  = $fields[22];
	    my $hg_intron	  = $fields[23];
	    my $rae_intron	  = $fields[24];
	    my $hg_5utr		  = $fields[25];
	    my $hg_3utr		  = $fields[26];
	    my @extra		  = split(",",$fields[27]);
	    my $notes		  = $fields[28];
	    
	    foreach my $syn (@synonyms){
		$entries{$syn}{GOI}=$goi;
		push @{$entries{$syn}{PATHWAY}}, @pathways;
		push @{$entries{$syn}{LITERATURE}}, @literature;
	    }
	    
	    $entries{$syn}{NAME}=$gene;
	    $entries{$syn}{TEX}="$goi\/$gio\.tex";
	    for (1..$igvs){
		push @{$entries{$syn}{IGV}},"$goi\/snapshots/$goi\_igv$_\.svg";
	    }
	    for (1..$extra[0]){
		push @{$entries{$syn}{UCSC}}"$goi\/snapshots/$goi\_ucsc$_\.eps";
	    }
	    $entries{$syn}{SASHIMI}="$goi\/snapshots/$goi\_sashimi.svg" if ($sashimi > 0);
	    $entries{$syn}{NOTES}=$notes;
	    $entries{$syn}{EXTRA}=$extra;
	    $entries{$syn}{CUFFLINKS}="Filled by Martin";
	    $entries{$syn}{READS}="Filled by Markus";
	}
    }
#    elsif ($filetoparse =~ /.xlsx/){
#	print STDERR "Parsing Excel sheet $filetoparse!\n";
#	my $parser = Spreadsheet::ParseXLSX->new;
#	my $workbook = $parser->parse("$filetoparse");
#	print STDERR $workbook->get_filename()."\n";
#	
#	if ( !defined $workbook ) {
#	    die $parser->error(), ".\n";
#	}
#	for my $worksheet ( $workbook->worksheets() ) {
#	    next unless ($worksheet->get_name() eq 'APG' || $worksheet->getname() eq 'GOI');
#	    my ( $row_min, $row_max ) = $worksheet->row_range();
#	    my ( $col_min, $col_max ) = $worksheet->col_range();
#
#	    for my $row ( $row_min .. $row_max ) {
#		for my $col ( $col_min .. $col_max ) {
#
#		    my $cell = $worksheet->get_cell( $row, $col );
#		    next unless $cell;
#		    print "Row, Col    = ($row, $col)\n";
#		    print "Value       = ", $cell->value(),       "\n";
#		    print "Unformatted = ", $cell->unformatted(), "\n";
#		    print "\n";
#		}
#	    }
#	}
#
#
#	my $excel = Spreadsheet::XLSX -> new ($file);
#	foreach my $sheet (@{$excel -> {Worksheet}}) {	    
#	    next unless ($sheet->{Name} eq 'APG' || $sheet->{Name} eq 'GOI');
#	    printf("Sheet: %s\n", $sheet->{Name});    
#	    $sheet -> {MaxRow} ||= $sheet -> {MinRow};	    
#	    foreach my $row ( $sheet -> {MinRow} .. $sheet -> {MaxRow}) {
#		print STDERR "ROW: $row\n";
#               $sheet -> {MaxCol} ||= $sheet -> {MinCol};                
#                foreach my $col ($sheet -> {MinCol} ..  $sheet -> {MaxCol}) {		    
#		    print STDERR "COL: $col\n";
#		    my $cell = $sheet -> {Cells} [$row] [$col];		    
#		    if ($cell) {
#			printf("( %s , %s ) => %s\n", $row, $col, $cell -> {Val});
#		    }	    
#                }		
#	    }	    
#	}
    }
    chdir ($odir) or die $!;
}

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
