#!/usr/bin/perl
### To use this Script replace all semicolons in fields with commas,
### then all line breaks in fields with nothing,
### then save as semicolon separated list and have fun parsing
### 
### Script supplementer.pl;
### Last changed Time-stamp: <2014-11-29 22:55:19 fall> by joerg

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
my $wdir = cwd();
my $today = DateTime->now->strftime('%d%m%Y');
my @csvs = ($goil,$apgl);
push @csvs, split(",",$peakl);
$odir =	 "Supplements_$today/" unless $odir;
$dir  =~ s/ //g;
$odir =~ s/ //g;
($dir) or die "No working directory chosen!\n";

my $pid = $$;
(my $job = `cat /proc/$pid/cmdline`)=~ s/\0/ /g;
print STDERR $job,"\n";

###############
###Main Stuff
###############

if (!-d $odir){
    print STDERR "Creating output directory $odir!\n";
    make_path ($odir) or die "Error creating directory: $odir";
}

my %genes; 
foreach my $file (@csvs){
    print STDERR "Reading input from $file!\n";
#### Read fields for html from file
#    my ($wdir,$filetoparse) = split("\/",$file,2);
    chdir ($dir) or die $!;
    if ($file =~ /goi|apg/i){
	%genes = %{read_tables($file,\%genes)};
    }
    else{
	%genes = %{parse_expression($file,\%genes)};
    }
    chdir ($wdir);
}

my ($html_destination_path) = ($odir);
make_supplements(\%genes,$html_destination_path);

sub make_supplements{
    my %gois = %{$_[0]};
#    print Dumper (\%gois);
    #check arguments
    die ("ERROR $html_destination_path does not exist\n") unless (-d $html_destination_path);
#    die ("ERROR no URL (network location) provided") unless(defined $base_URL);
    chdir($odir) or die "$!";

    my $template_path = "supplementer/scripts/template";

    #template definition
    my $template = Template->new({
	INCLUDE_PATH => ["$template_path"],
	RELATIVE=>1,
				 });
    #ensure that base_URL ends with slash
    #$base_URL =~ s!/*$!/!;  
    my $logPath="$html_destination_path/Log";
    
    #create html directory structure
    my @genelist;
    my @parseit = ('GOI', 'APG', 'EXPRESSION');
    foreach my $gene( keys %gois ){
	push @genelist, $gene;
	foreach my $from (@parseit){
#	    print STDERR "$from\t$gene\n";
	    next unless (defined $gois{$gene}{$from}{ID});
	    my $goi = $gois{$gene}{$from}{ID};	
	    #construct gene of interest goi.html
	    my $goi_path = $odir . "/goi.html";
	    my $goi_file = $goi.".html";
	    my $name = $gene;
	    my ($cufflinks, $maxy);
	    foreach my $sample (keys %{$gois{$gene}{$from}{CUFFLINKS}} ){
		$cufflinks .= $sample.":";
		$cufflinks .= join(",",@{$gois{$gene}{$from}{CUFFLINKS}{$sample}}) if (defined $gois{$gene}{$from}{CUFFLINKS});
		$maxy .= $sample.":";
		$maxy .= join(",",@{$gois{$gene}{$from}{PEAKS}{$sample}}) if (defined $gois{$gene}{$from}{PEAKS});
	    }
	    ($cufflinks, $maxy) = ('NA','NA') unless defined ($cufflinks && $maxy);
            my $igv = image_entry(${$gois{$gene}{$from}{IGV}}[0],$dir,$odir);
            my $sashimi = image_entry(${$gois{$gene}{$from}{SASHIMI}}[0],$dir,$odir);
            my $ucsc = image_entry(${$gois{$gene}{$from}{UCSC}}[0],$dir,$odir);
	    my $goi_vars = 
	    {   
		name => $gois{$gene}{$from}{NAME},
		synonyms => join(",",$gois{$gene}{$from}{SYNONYMS}),
		goiid => $gois{$gene}{$from}{ID},
		textxt => $gois{$gene}{$from}{NOTES},
		igv => $igv,
		sashimi => $sashimi,
		ucsc => $ucsc,
		additionalplots => $gois{$gene}{$from}{EXTRA},
		cufflinks => $cufflinks,
		maxy => $maxy 
	    };
	    print Dumper(\$goi_vars); 
	    $template->process($goi_file,$goi_vars,$goi_path) || die "Template process failed: ", $template->error(), "\n";	
	}
	
	#construct index.hmtl
#	my $index_path = $html_destination_path. "/index.html";
#	my $index_file = 'index.html';
#	my $index_vars = index_entry();
#	$template->process($index_file,$index_vars,$index_path) || die "Template process failed: ", $template->error(), "\n";
    }
    chdir($wdir) or die "$!";
}

sub parse_expression{
#HG19 Mock vs Ebola 
#HG19 Mock vs MARV
#HG19 EBOV vs MARV
#RAE Mock vs EBOV
#RAE Mock vs MARV
#RAE EBOV vs MARV
#HG19 EBOV vs RAE EBOV
#HG19 MARV vs RAE MARV
    my $filetoparse = $_[0];
    my $sample = (split(/\./,$filetoparse))[0];
    my %entries	    = %{$_[1]};
    print STDERR "Expression parsing $sample!\n";
    open (LIST,"<:gzip(autopop)","$filetoparse");
    while(<LIST>){
	next if($_ =~ /^#/);
	my $line  = $_;
	my ($gene, $mb3, $mb7, $mb23, $eb3, $eb7, $eb23, $l3, $l7, $l23, $max, $mp3, $mp7, $mp23, $ep3, $ep7, $ep23) = split(/\t/,$line);
#	print STDERR $gene,"\n";
	my $goto;
	if (defined $entries{$gene}{GOI}){
	    $goto = "GOI";
	}
	elsif(defined $entries{$gene}{APG}){
	    $goto = "APG";
	}
	else{
	    $goto = "EXPRESSION";
	}
	push @{$entries{$gene}{$goto}{CUFFLINKS}{$sample}}, ($mb3, $mb7, $mb23, $eb3, $eb7, $eb23);
	push @{$entries{$gene}{$goto}{LOGEXPRESSION}{$sample}}, ($l3, $l7, $l23);
	push @{$entries{$gene}{$goto}{MAX}{$sample}}, $max;
	push @{$entries{$gene}{$goto}{PEAKS}{$sample}}, ($mp3, $mp7, $mp23, $ep3, $ep7, $ep23);
    }
    return (\%entries);
}

sub read_tables{
    my $filetoparse = $_[0];
    my %entries	    = %{$_[1]};
    my @again;
    print STDERR "Parsing $filetoparse!\n";
    open (my $list,"<:gzip(autopop)","$filetoparse");
    my @process = <$list>;

    if ($filetoparse =~ /goi/i){
	print STDERR "Processing GIO List!\n";
	foreach(@process){
	    next unless ($_ =~ /.goi./);
	    (my $line	  = $_) =~ s/,w+//g;
	    my @fields		  = split(/\;/,$line);
	    my $goi		  = $fields[0];
	    my $hacker		  = $fields[1];
	    my $gene		  = $fields[2];
	    next if (defined $entries{$gene}{GOI}{ID});

	    my $duplicate	  = $fields[3];   
	    my @synonyms	  = split(/[,\s]+/,$fields[4]);
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
		
		$entries{$gene}{GOI}{ID}  =	$goi;
		$entries{$gene}{GOI}{SYNONYMS} = @synonyms;
		push @{$entries{$gene}{GOI}{PATHWAY}}, @pathways;
		push @{$entries{$gene}{GOI}{LITERATURE}}, @literature;   
		$entries{$gene}{GOI}{NAME} =	$gene;
		$entries{$gene}{GOI}{TEX}  =	"$goi\/$goi\.tex";
		if ($igvs == 1){
		    push @{$entries{$gene}{GOI}{IGV}},"$goi\/snapshots/$goi\_igv.svg";
		}
		elsif ($igvs == 0){
		    push @{$entries{$gene}{GOI}{IGV}},"NONE";
		}
		else{
		    for (1..$igvs){
			push @{$entries{$gene}{GOI}{IGV}},"$goi\/snapshots/$goi\_igv$_\.svg";
		    }
		}
		if ($ucscs == 1){
		    push @{$entries{$gene}{GOI}{UCSC}},"$goi\/snapshots/$goi\_ucsc.svg";
		}
		elsif ($ucscs == 0){
		    push @{$entries{$gene}{GOI}{UCSC}},"NONE";
		}
		else{
		    for (1..$ucscs){
			push @{$entries{$gene}{GOI}{UCSC}},"$goi\/snapshots/$goi\_ucsc$_\.svg";
		    }
		}
		if ($sashimi == 1){
		    push @{$entries{$gene}{GOI}{SASHIMI}},"$goi\/snapshots/$goi\_sashimi.svg";
		}
		elsif ($sashimi == 0){
		    push @{$entries{$gene}{GOI}{SASHIMI}},"NONE";
		}
		else{
		    for (1..$sashimi){
			push @{$entries{$gene}{GOI}{SASHIMI}},"$goi\/snapshots/$goi\_sashimi$_\.svg";
		    }
		}
		$entries{$gene}{GOI}{NOTES}	   = $notes || 'NA';
		$entries{$gene}{GOI}{EXTRA}	   = $extra || 'NA';
		
		foreach my $syn (@synonyms){
		    next if ($syn eq $gene);
		    $entries{$syn}{GOI} = $entries{$gene}{GOI};		    
		}   
	    }
	    else{
		foreach my $syn (@synonyms){
		    next if ($syn eq $gene);
#		    print STDERR "GOI: $gene is a Duplicate of $syn but this has no been processed yet\n" unless ($entries{$syn});
		    push @process, $line unless ($entries{$syn});
		    $entries{$gene}{GOI} = $entries{$syn}{GOI} if ($entries{$syn}{GOI});
		}
	    }
	}
    }
    elsif ($filetoparse =~ /apg/){
	print STDERR "Processing APG List!\n";
	foreach(@process){
	    next unless ($_ =~ /.apg./);
	    (my $line	  = $_) =~ s/,w+//g;
	    my @fields		  = split(/\;/,$line);
	    my $apg		  = $fields[0];
	    my $hacker		  = $fields[1];
	    my $gene		  = $fields[2];
	    next if (defined $entries{$gene}{APG}{ID});

	    my $duplicate	  = $fields[3];   
	    my @synonyms	  = split(/[,\s]+/,$fields[4]);
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
		
		$entries{$gene}{APG}{ID}  = $apg;
		$entries{$gene}{APG}{SYNONYMS} = @synonyms;
		push @{$entries{$gene}{APG}{PATHWAY}}, @pathways;
		push @{$entries{$gene}{APG}{LITERATURE}}, @literature;   
		$entries{$gene}{APG}{NAME} =	$gene;
		$entries{$gene}{APG}{TEX}  =	"$apg\/$apg\.tex";
		if ($igvs == 1){
		    push @{$entries{$gene}{APG}{IGV}},"$apg\/snapshots/$apg\_igv.svg";
		}
		elsif ($igvs == 0){
		    push @{$entries{$gene}{APG}{IGV}},"NONE";
		}
		else{
		    for (1..$igvs){
			push @{$entries{$gene}{APG}{IGV}},"$apg\/snapshots/$apg\_igv$_\.svg";
		    }
		}
		if ($ucscs == 1){
		    push @{$entries{$gene}{APG}{UCSC}},"$apg\/snapshots/$apg\_ucsc.svg";
		}
		elsif ($ucscs == 0){
		    push @{$entries{$gene}{APG}{UCSC}},"NONE";
		}
		else{
		    for (1..$ucscs){
			push @{$entries{$gene}{APG}{UCSC}},"$apg\/snapshots/$apg\_ucsc$_\.svg";
		    }
		}
		if ($sashimi == 1){
		    push @{$entries{$gene}{APG}{SASHIMI}},"$apg\/snapshots/$apg\_sashimi.svg";
		}
		elsif ($sashimi == 0){
		    push @{$entries{$gene}{APG}{SASHIMI}},"NONE";
		}
		else{
		    for (1..$sashimi){
			push @{$entries{$gene}{APG}{SASHIMI}},"$apg\/snapshots/$apg\_sashimi$_\.svg";
		    }
		}
		$entries{$gene}{APG}{NOTES}	   = $notes || 'NA';
		$entries{$gene}{APG}{EXTRA}	   = $extra || 'NA';
		
		foreach my $syn (@synonyms){
		    next if ($syn eq $gene);
		    $entries{$syn}{APG}=$entries{$gene}{APG};		    
		}   
	    }
	    else{
		foreach my $syn (@synonyms){
		    next if ($syn eq $gene);
#		    print STDERR "APG: $gene is a Duplicate of $syn but this has no been processed yet\n" unless ($entries{$syn});
		    push @process, $line unless ($entries{$syn});
		    $entries{$gene}{APG}=$entries{$syn}{APG} if ($entries{$syn}{APG});
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

sub index_entry{
    my $synonym = shift;
    my $goi = shift;
    my $goilink = $goi.".html";
    my $index_entry = "<tr><td><a href=\"$goilink\">$synonym</a></td></tr>";
    return $index_entry;
}

sub image_entry{
    my $file = shift;
    my $dir = shift;
    my $odir = shift;
    my @file = split ("/", $file);
    my $filename = $file[2];
    my $imagelink = $dir ."/". $file;
    my $thumblink = $odir ."thumbs/" . "$filename"; 
    `convert $imagelink -resize 150×150! $thumblink`;
    my $image_entry = "<a href=\"$imagelink\"><img src=\"$thumblink\"></a>";
    return $image_entry;
}
