#!/usr/bin/perl
### To use this Script replace all semicolons in fields with commas,
### then all line breaks in fields with nothing,
### then save as semicolon separated list and have fun parsing
### 
### Script supplementer.pl;
### Last changed Time-stamp: <2014-11-28 18:38:08 fall> by joerg

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
#use PerlIO::gzip;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use DateTime qw();
use Data::Dumper;

###############
###Variables
###############
my $VERBOSE = 0;
my ( $dir, $odir, $files, $outfiles, $goil, $apgl, $peakl );

###############
###Command Line Options
###############

pod2usage(-verbose => 0)
    unless GetOptions(
	"dir|d=s"    => \$dir,
	"odir|o=s"   => \$odir,
	"goi|g=s"    => \$goil,
	"apg|a=s"    => \$apgl,
	"peaks|p=s"  => \$peakl,
	"ofiles|t=s" => \$outfiles,
	"help|h"     => sub{pod2usage(-verbose => 1)},
	"man|m"	     => sub{pod2usage(-verbose => 2)},
	"verbose"    => sub{ $VERBOSE++ }
    );

$dir  =	 cwd() unless ($dir);
my $today = DateTime->now->strftime('%d%m%Y');
my @csvs = ($peakl,$goil,$apgl);
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

my %genes; 
foreach my $file (@csvs){
    print STDERR "Reading input!\n";
    die "File $file could not be openend!\n" unless (-e $file);
#### Read fields for html from file
    my ($wdir,$filetoparse) = split("\/",$file,2);
    chdir ($wdir) or die $!;
    
    if ($filetoparse =~ /goi|apg/i ){
	%genes = %{read_tables($filetoparse,\%genes)};
    }
    else{
	%genes = %{parse_expression($filetoparse,\%genes)};
    }
    chdir ($dir) or die $!;
}

chdir ($odir) or die $!;
make_supplements(\%genes,$html_destination_path,$base_URL,);

sub make_supplements{
    my %gois = %{$_[0]};
    #   foreach my $key( keys %gois ){
    #	print STDERR @{$gois{$key}{PEAKS}} if (defined $gois{$key}{PEAKS});
    #   }
    print Dumper (\%gois);
    #check arguments
    die ("ERROR $html_destination_path does not exist\n") unless (-d $html_destination_path);
    die ("ERROR no URL (network location) provided") unless(defined $base_URL);

    #ensure that base_URL ends with slash
    $base_URL =~ s!/*$!/!;  
    my $logPath="$html_destination_path/Log"
    #check program dependencies


    #create html directory structure
    
    #template definition
    my $template = Template->new({
	INCLUDE_PATH => ["$template_path"],
	RELATIVE=>1,
 				 });

    #construct index.hmtl
    my ($foo, $bar)= (0,0);
    my $index_path = $html_destination_path. "/index.html";
    my $index_file = 'index.html';
    my $index_vars =
    {
	genesofinterests => "$genesofinterests"
        
    };
    $template->process($index_file,$index_vars,$index_path) || die "Template process failed: ", $template->error(), "\n";

    #construct gene of interest goi.html
    my $goi_path = $assembly_hub_directory . "/goi.html";
    my $goi_file = 'goi.html';
    my $goi_vars =
    {
	name => "$name",
        synonyms => "$synonyms",
       	goiid => "$goiid",
       	textxt => "$textxt",
       	igv => "$igv",
        sashimi => "$sashimi",
        ucsc => "$ucsc",
        additionalplots => "$additionalplots",
        cufflinks => "$cufflinks",
        maxy => "$maxy"    
    };
    $template->process($goi_file,$goi_vars,$goi_path) || die "Template process failed: ", $template->error(), "\n";
}

sub parse_expression{
    my $filetoparse = $_[0];
    my %entries	    = %{$_[1]};
    print STDERR "Expression parsing $filetoparse!\n";
    open (LIST,"<:gzip(autopop)","$filetoparse");
    while(<LIST>){
	next if($_ =~ /^#/);
	my $line  = $_;
	my ($gene, $mb3, $mb7, $mb23, $eb3, $eb7, $eb23, $l3, $l7, $l23, $max, $mp3, $mp7, $mp23, $ep3, $ep7, $ep23) = split(/\t/,$line);
#	print STDERR $gene,"\n";
	push @{$entries{$gene}{CUFFLINKS}}, ($mb3, $mb7, $mb23, $eb3, $eb7, $eb23);
	push @{$entries{$gene}{LOGEXPRESSION}}, ($l3, $l7, $l23);
	push @{$entries{$gene}{MAX}}, $max;
	push @{$entries{$gene}{PEAKS}}, ($mp3, $mp7, $mp23, $ep3, $ep7, $ep23);
    }
    return (\%entries);
}

sub read_tables{
    my $filetoparse = $_[0];
    my %entries	    = %{$_[1]};
    print STDERR "Parsing $filetoparse!\n";
    open (LIST,"<:gzip(autopop)","$filetoparse");
    if ($filetoparse =~ /goi/i){
	print STDERR "Processing GIO List!\n";
	while(<LIST>){
	    next unless ($_ =~ /.goi./);
	    (my $line	  = $_) =~ s/,w+//g;
	    my @fields		  = split(/\;/,$line);
	    my $goi		  = $fields[0];
	    my $hacker		  = $fields[1];
	    my $gene		  = $fields[2];
	    next if (defined $entries{$gene}{GOI});

	    my $duplicate	  = $fields[3];   
	    my @synonyms	  = split(",",$fields[4]);
	    push @synonyms, $gene unless ($synonyms[0]);
	    if ($duplicate eq '' || $duplicate == 0 || ($duplicate && $fields[8] ne '')){
		my @pathways	      = split(",",$fields[5]) if ($fields[5] ne '');
		push @pathways, 'Unknown' unless ($pathways[0]);
		my @literature	      = split(",",$fields[6]) if ($fields[6] ne '');
		push @literature, 'Unknown' unless ($literature[0]);
		my $igvs	      = (split(/[:,\s\/]+/,$fields[7],2))[0] || '0';
		$igvs = 0 if ($igvs =~ /todo/i);
		my $mock	      = $fields[8];
		my $ebola	      = $fields[9];
		my $marburg	      = $fields[10];
		my $profile	      = $fields[11];
		my $homolog	      = $fields[12];
		my $raemock	      = $fields[13];
		my $raeebola	      = $fields[14];
		my $raemarburg	      = $fields[15];
		my $raeprofile	      = $fields[16];
		my $ucscs	      = (split(/[:,\s\/]+/,$fields[17],2))[0] || '0';
		$ucscs = 0 if ($ucscs =~ /todo/i);
		my $ucsc_conservation = $fields[18];
		my $hg_SNPs	      = $fields[19];
		my $rae_SNPs	      = $fields[20];
		my $sashimi	      = (split(/[:,\s\/]+/,$fields[21],2))[0] || '0';
		$sashimi = 0 if ($sashimi =~ /todo/i);
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
		
		foreach my $syn (@synonyms){
		    next if ($syn eq $gene);
		    $entries{$syn} = $entries{$gene};		    
		}   
	    }
	    else{
		foreach my $syn (@synonyms){
		    next if ($syn eq $gene);
		    print STDERR "GOI: $gene is a Duplicate of $syn but this has no been processed yet\n" unless ($entries{$syn});
		    $entries{$gene} = $entries{$syn} if ($entries{$syn});
		}
	    }
	}
    }
    elsif ($filetoparse =~ /apg/){
	print STDERR "Processing APG List!\n";
	while(<LIST>){
	    next unless ($_ =~ /.apg./);
	    (my $line	  = $_) =~ s/,w+//g;
	    my @fields		  = split(/\;/,$line);
	    my $apg		  = $fields[0];
	    my $hacker		  = $fields[1];
	    my $gene		  = $fields[2];
	    next if (defined $entries{$gene}{APG});

	    my $duplicate	  = $fields[3];   
	    my @synonyms	  = split(",",$fields[4]);
	    push @synonyms, $gene unless ($synonyms[0]);
	    if ($duplicate eq '' || $duplicate == 0){
		my @pathways	      = split(/[:,]+/,$fields[5]) if ($fields[5] ne '');
		push @pathways, 'Unknown' unless ($pathways[0]);
		my @literature	      = split(/[:,]+/,$fields[6]) if ($fields[6] ne '');
		push @literature, 'Unknown' unless ($literature[0]);
		my $igvs	      = (split(/[:,\s\/]+/,$fields[7],2))[0] || '0';
		$igvs = 0 if ($igvs =~ /todo/i);
		my $scale             = $fields[8] || '0';
		my $mock	      = $fields[9];
		my $ebola	      = $fields[10];
		my $marburg	      = $fields[11];
		my $profile	      = $fields[12];
		my $homolog	      = $fields[13];
		my $raemock	      = $fields[14];
		my $raeebola	      = $fields[15];
		my $raemarburg	      = $fields[16];
		my $raeprofile	      = $fields[17];
		my $ucscs	      = (split(/[:,\s\/]+/,$fields[17],2))[0] || '0';
		$ucscs = 0 if ($ucscs =~ /todo/i);
		my $ucsc_conservation = $fields[19];
		my $hg_SNPs	      = $fields[20];
		my $rae_SNPs	      = $fields[21];
		my $sashimi	      = (split(/[:,\s\/]+/,$fields[21],2))[0] || '0';
		$sashimi = 0 if ($sashimi =~ /todo/i);
		my $hg_intron	      = $fields[23];
		my $rae_intron	      = $fields[24];
		my $hg_5utr	      = $fields[25];
		my $hg_3utr	      = $fields[26];
		my $extra	      = $fields[27];
		my $notes	      = $fields[28];
		
		$entries{$gene}{APG}  =	$apg;
		push @{$entries{$gene}{PATHWAY}}, @pathways;
		push @{$entries{$gene}{LITERATURE}}, @literature;   
		$entries{$gene}{NAME} =	$gene;
		$entries{$gene}{TEX}  =	"$apg\/$apg\.tex";
		if ($igvs == 1){
		    push @{$entries{$gene}{IGV}},"$apg\/snapshots/$apg\_igv.svg";
		}
		elsif ($igvs == 0){
		    push @{$entries{$gene}{IGV}},"NONE";
		}
		else{
		    for (1..$igvs){
			push @{$entries{$gene}{IGV}},"$apg\/snapshots/$apg\_igv$_\.svg";
		    }
		}
		if ($ucscs == 1){
		    push @{$entries{$gene}{UCSC}},"$apg\/snapshots/$apg\_ucsc.svg";
		}
		elsif ($ucscs == 0){
		    push @{$entries{$gene}{UCSC}},"NONE";
		}
		else{
		    for (1..$ucscs){
			push @{$entries{$gene}{UCSC}},"$apg\/snapshots/$apg\_ucsc$_\.svg";
		    }
		}
		$entries{$gene}{SASHIMI}   = "$apg\/snapshots/$apg\_sashimi.svg" if ($sashimi > 0);
		$entries{$gene}{NOTES}	   = $notes;
		$entries{$gene}{EXTRA}	   = $extra;
		
		foreach my $syn (@synonyms){
		    next if ($syn eq $gene);
		    $entries{$syn}=$entries{$gene};		    
		}   
	    }
	    else{
		foreach my $syn (@synonyms){
		    next if ($syn eq $gene);
		    print STDERR "APG: $gene is a Duplicate of $syn but this has no been processed yet\n" unless ($entries{$syn});
		    $entries{$gene}=$entries{$syn} if ($entries{$syn});
		}
	    }
	}	
    }
    return (\%entries);
}

sub unique_array{

    my $arrayref = shift;
    my @array = @{$arrayref};

    my %unique = ();
    foreach my $item (@array)
    {
        $unique{$item} ++;
    }
    my @arrayuid = sort {$a cmp $b} keys %unique;

    return(\@arrayuid);
}
