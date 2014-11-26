#!/usr/bin/perl

### To use this Script replace all semicolons in fields with commas,
### then all line breaks in fields with nothing,
### then save as semicolon separated list and have fun parsing
### 
### Script supplementer.pl;
### Last changed Time-stamp: <2014-11-26 17:21:58 fall> by joerg
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
use Data::Dumper;

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

    if ($filetoparse =~ /[.csv|.ssv|.tsv]/ ){
#hg19.goi.00001;manja;MD2;0;LY96;;;1;00;00;00;0;0-;00;00;00;0;1;11110;0;0;1;0;0;0;0;0;0;#VALUE!;1;1;0;0;1;.;CHECKED
	print STDERR "Parsing $filetoparse!\n";
	open (LIST,"<:gzip(autopop)","$filetoparse");
	while(<LIST>){
	    next unless ($_ =~ /.goi./ || $_ =~ /.apk./);
#	    print STDERR $_;
	    (my $line	  = $_) =~ s/,w+//g;
	    my @fields		  = split(/\;/,$line);
#	    print STDERR @fields,"\n";
	    my $goi		  = $fields[0];
	    my $hacker		  = $fields[1];
	    my $gene		  = $fields[2];
	    my $duplicate	  = $fields[3];   
	    my @synonyms	  = split(",",$fields[4]);
	    push @synonyms, $gene unless ($synonyms[0]);
	    if ($duplicate eq '' || $duplicate == 0){
		my @pathways	      = split(",",$fields[5]) if ($fields[5] ne '');
		push @pathways, 'Unknown' unless ($pathways[0]);
		my @literature	      = split(",",$fields[6]) if ($fields[6] ne '');
		push @literature, 'Unknown' unless ($literature[0]);
		my $igvs	      = $fields[7] || '0';
		my $mock	      = $fields[8];
		my $ebola	      = $fields[9];
		my $marburg	      = $fields[10];
		my $profile	      = $fields[11];
		my $homolog	      = $fields[12];
		my $raemock	      = $fields[13];
		my $raeebola	      = $fields[14];
		my $raemarburg	      = $fields[15];
		my $raeprofile	      = $fields[16];
		my $ucscs	      = $fields[17] || '0';
		my $ucsc_conservation = $fields[18];
		my $hg_SNPs	      = $fields[19];
		my $rae_SNPs	      = $fields[20];
		my $sashimi	      = $fields[21] || '0';
		my $hg_intron	      = $fields[22];
		my $rae_intron	      = $fields[23];
		my $hg_5utr	      = $fields[24];
		my $hg_3utr	      = $fields[25];
		my $extra	      = $fields[26];
		my $notes	      = $fields[27];
		
		$entries{$gene}{GOI}  =	$goi;
		push @{$entries{$gene}{PATHWAY}}, @pathways;
		push @{$entries{$gene}{LITERATURE}}, @literature;   
		$entries{$gene}{NAME} =	$gene;
		$entries{$gene}{TEX}  =	"$goi\/$goi\.tex";
		if ($igvs == 1){
		    push @{$entries{$gene}{IGV}},"$goi\/snapshots/$goi\_igv.svg";
		}
		elsif ($igvs == 0){
		    push @{$entries{$gene}{IGV}},"NONE";
		}
		else{
		    for (1..$igvs){
			push @{$entries{$gene}{IGV}},"$goi\/snapshots/$goi\_igv$_\.svg";
		    }
		}
		if ($ucscs == 1){
		    push @{$entries{$gene}{UCSC}},"$goi\/snapshots/$goi\_ucsc.svg";
		}
		elsif ($ucscs == 0){
		    push @{$entries{$gene}{UCSC}},"NONE";
		}
		else{
		    for (1..$ucscs){
			push @{$entries{$gene}{UCSC}},"$goi\/snapshots/$goi\_ucsc$_\.svg";
		    }
		}
		$entries{$gene}{SASHIMI}   = "$goi\/snapshots/$goi\_sashimi.svg" if ($sashimi > 0);
		$entries{$gene}{NOTES}	   = $notes;
		$entries{$gene}{EXTRA}	   = $extra;
		$entries{$gene}{CUFFLINKS} = "Filled by Martin";
		$entries{$gene}{READS}	   = "Filled by Markus";
		
		foreach my $syn (@synonyms){
		    next if ($syn eq $gene);
		    $entries{$syn}=$entries{$gene};		    
		}   
	    }
	    else{
		foreach my $syn (@synonyms){
		    $entries{$gene}=$entries{$syn} if ($entries{$syn});
		}
	    }
	}
    }
    
    chdir ($odir) or die $!;
    make _supplements(\%entries);
}

sub make_supplements{
    my %gois = %{$_};
    
    #check arguments
    die ("ERROR $html_destination_path does not exist\n") unless (-d $assembly_hub_destination_path);
    die ("ERROR no URL (network location) provided") unless(defined $base_URL);
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
