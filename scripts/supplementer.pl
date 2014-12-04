#!/usr/bin/perl
### To use this Script replace all semicolons in fields with commas,
### then all line breaks in fields with nothing,
### then save as semicolon separated list and have fun parsing
### 
### Script supplementer.pl;
### Last changed Time-stamp: <2014-12-04 11:59:33 fall> by joerg

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

if (!-d "$odir/thumbs"){
    print STDERR "Creating output directory $odir/thumbs!\n";
    make_path ("$odir/thumbs") or die "Error creating directory: $odir/thumbs $!";
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
    elsif($file =~ /mock|ebov|marv/i){
	%genes = %{parse_expression($file,\%genes)};
    }
    elsif($file =~ /timepoints/i){
	%genes = %{parse_time($file,\%genes)};
    }
    elsif($file =~ /hg19/i && $file =~ /rae/i){
	%genes = %{parse_comparison($file,\%genes)};
    }
    chdir ($wdir);
}

my $html_destination_path = join("/",$wdir,$odir);
make_supplements(\%genes,$html_destination_path);

sub make_supplements{
    my %gois = %{$_[0]};
#    print Dumper(\%gois);
    #check arguments
    die ("ERROR $html_destination_path does not exist\n") unless (-d $html_destination_path);
#    die ("ERROR no URL (network location) provided") unless(defined $base_URL);
    chdir($odir) or die "$!";
    my $template_path = "$wdir/supplementer/scripts/template";

    #template definition
    my $template = Template->new({
	INCLUDE_PATH => ["$template_path"],
	RELATIVE=>0
				 });
    #ensure that base_URL ends with slash
    #$base_URL =~ s!/*$!/!;  
    my $logPath="$html_destination_path/Log";
    
    #create html directory structure
    my %genelist;
    my @parseit = ('GOI', 'APG', 'EXPRESSION');
#    my @sorted_genes = 	map { $_->[0] } sort { $a->[1] cmp $b->[1] } map { [ $_, uc($_) ] } keys %gois;
    foreach my $gene( sort {lc($a) cmp lc($b)} keys %gois ){
#    foreach my $gene( @sorted_genes ){
	foreach my $from (@parseit){
	    next unless (defined $gois{$gene}{$from}{ID});
	    my $goi = $gois{$gene}{$from}{ID};	
	    #construct gene of interest goi.html
	    my $goi_path = "$goi.html";
	    my $goi_file = "goi.html";
	    my $name = $gene;
	    $genelist{$name}=$goi;
### Parse Expression
	    my (@samp, @condi, @deg, @max, @fold) = ();
	    foreach my $sample (keys %{$gois{$gene}{$from}{CUFFLINKS}} ){
		my $fold_change = join(",",@{$gois{$gene}{$goto}{LOGEXPRESSION}{$sample}{LOG}}) if (defined $gois{$gene}{$goto}{LOGEXPRESSION}{$sample}{LOG});
		push @fold, $fold_change;
		push @samp, $sample;
		foreach my $condition (keys %{$gois{$gene}{$from}{CUFFLINKS}{$sample}} ){
		    my $cufflinks = join(",",@{$gois{$gene}{$from}{CUFFLINKS}{$sample}{$condition}}) if (defined $gois{$gene}{$from}{CUFFLINKS}{$sample}{$condition});
		    my $maxy = join(",",@{$gois{$gene}{$from}{PEAKS}{$sample}{$condition}}) if (defined $gois{$gene}{$from}{PEAKS}{$sample}{$condition});
		    push @deg, $cufflinks;
		    push @max, $maxy;
		    push @condi, $condition;
		}
	    }

### Parse Comparison
	    my (@csamp, @ccondi, @cdeg, @cmax, @cfold) = ();
	    foreach my $sample (keys %{$gois{$gene}{$from}{HB}} ){
		my $cfold_change = join(",",@{$entries{$gene}{$goto}{CLOGEXPRESSION}{$sample}{LOG}});
		push @cfold, $cfold_change;
		push @csamp, $sample;
		foreach my $condition (keys %{$gois{$gene}{$from}{HB}{$sample}} ){
		    my $cufflinks = join(",",@{$gois{$gene}{$from}{HB}{$sample}{$condition}}) if (defined $gois{$gene}{$from}{HB}{$sample}{$condition});
		    push @cdeg, $cufflinks;	    
		    push @ccondi, $condition;
		}
		foreach my $condition (keys %{$gois{$gene}{$from}{CPEAKS}{$sample}} ){
		    my $maxy = join(",",@{$gois{$gene}{$from}{CPEAKS}{$sample}{$condition}}) if (defined $gois{$gene}{$from}{CPEAKS}{$sample}{$condition});
		    push @cmax, $maxy;
		}
	    }

### Parse Timepoints
	    my (@tsamp, @tcondi, @tdeg, @tmax, @tfold) = ();
	    foreach my $sample (keys %{$gois{$gene}{$from}{MEV}} ){
		push @tsamp, $sample;
		foreach my $condition (keys %{$gois{$gene}{$from}{MEV}{$sample}} ){
		    my $cufflinks = join(",",@{$gois{$gene}{$from}{MEV}{$sample}{$condition}}) if (defined $gois{$gene}{$from}{MEV}{$sample}{$condition});
		    my $maxy = join(",",@{$gois{$gene}{$from}{TPEAKS}{$sample}{$condition}}) if (defined $gois{$gene}{$from}{PEAKS}{$sample}{$condition});
		    my $fold_change = join(",",@{$gois{$gene}{$goto}{TLOGEXPRESSION}{$sample}{LOG}}) if (defined $gois{$gene}{$goto}{TLOGEXPRESSION}{$sample}{LOG});
		    push @tdeg, $cufflinks;
		    push @tmax, $maxy;
		    push @tfold, $fold_change;
		    push @tcondi, $condition;
		}
	    }
	    
            my $igv = image_entry(${$gois{$gene}{$from}{IGV}}[0],$dir,$odir);
            my $sashimi = image_entry(${$gois{$gene}{$from}{SASHIMI}}[0],$dir,$odir);
            my $ucsc = image_entry(${$gois{$gene}{$from}{UCSC}}[0],$dir,$odir);
            my $tex_link = link_entry($gois{$gene}{$from}{TEX},$dir);
	    (my $syn = join(",",@{$gois{$gene}{$from}{SYNONYMS}})) =~ s/, /,/g;
	    my $goi_vars = 
	    {   
		name		  => $gois{$gene}{$from}{NAME},
		synonyms	  => $syn,
		goiid		  => $gois{$gene}{$from}{ID},
		textxt		  => $tex_link,
		igv		  => $igv,
		sashimi		  => $sashimi,
		ucsc		  => $ucsc,
		additionalplots	  => $gois{$gene}{$from}{EXTRA},
##sample1
		sample		  => $samp[0],
##sample 1 condition1
		condone		  => $condi[0],
		maxy		  => $maxy[0],
		cufflinks	  => $deg[0],
##sample1 condition2
		condtwo		  => $condi[1],
 		maxy_two	  => $maxy[1],
		cufflinks_two	  => $deg[1],
##sample2
		sampleone	  => $samp[1],
##sample 2 condition1
		condoone	  => $condi[2],
		maxyo		  => $maxy[2],
		cufflinkso	  => $deg[2],
##sample 2 condition2
		condotwo	  => $condi[3],
		maxyo_two	  => $maxy[3],
		cufflinkso_two	  => $deg[3],
##sample 3
		sampletwo	  => $samp[2],
##sample 3 condition1
		condtone	  => $condi[4],
		maxyt		  => $maxy[4],
		cufflinkst	  => $deg[4],
##sample 3 condition2
		condttwo	  => $condi[5],
		maxyo_two	  => $maxy[5],
		cufflinkst_two	  => $deg[5],
###Comparison
##sample1
		csample		  => $csamp[0],
##sample 1 condition1
		ccondone	  => $ccondi[0],
		cmaxy		  => $cmaxy[0],
		ccufflinks	  => $cdeg[0],
##sample1 condition2
		ccondtwo	  => $ccondi[1],
 		cmaxy_two	  => $cmaxy[1],
		ccufflinks_two	  => $cdeg[1],
##sample2
		csampleone	  => $csample[1],
##sample 2 condition1
		ccondoone	  => $ccondi[2],
		cmaxyo		  => $cmaxy[2],
		ccufflinkso	  => $cdeg[2],
##sample 2 condition2
		ccondotwo	  => $ccondi[3],
		cmaxyo_two	  => $cmaxy[3],
		ccufflinkso_two	  => $cdeg[3],
### Timepoints
##sample1
		tsample		  => $tsample[0],
##sample 1 condition1
		tcondone	  => $tcondi[0],
		tmaxy		  => $tmaxy[0],
		tcufflinks	  => $tdeg[0],
##sample1 condition2
		tcondtwo	  => $tcondi[1],
 		tmaxy_two	  => $tmaxy[1],
		tcufflinks_two	  => $tdeg[1],
##sample1 condition3
		tcondthree	  => $tcondi[2],
 		tmaxy_three	  => $tmaxy[2],
		tcufflinks_three  => $tdeg[2],
##sample2
		tsampleone	  => $tsample[1],
##sample 2 condition1
		tcondoone	  => $tcondi[3],
		tmaxyo		  => $tmaxy[3],
		tcufflinkso	  => $tdeg[3],
##sample 2 condition2
		tcondotwo	  => $tcondi[4],
		tmaxyo_two	  => $tmaxy[4],
		tcufflinkso_two	  => $tdeg[4],
##sample 2 condition3
		tcondthree	  => $tcondi[5],
		tmaxyo_three	  => $tmaxy[5],
		tcufflinkso_three => $tdeg[5],
##sample 3
		tsampletwo	  => $tsample[2],
##sample 3 condition1
		tcondtone	  => $tcondi[6],
		tmaxyt		  => $tmaxy[6],
		tcufflinkst	  => $tdeg[6],
##sample 3 condition2
		tcondttwo	  => $tcondi[7],
		tmaxyt_two	  => $tmaxy[7],
		tcufflinkst_two	  => $tdeg[7],		
##sample 3 condition3
		tcondtthree	  => $tcondi[8],
		tmaxyt_three	  => $tmaxy[8],
		tcufflinkst_three => $tdeg[8]
	    };
	    $template->process($goi_file,$goi_vars,$goi_path) || die "Template process failed: ", $template->error(), "\n";	
	}
    }
    #construct index.hmtl
    my $index_path = $html_destination_path. "/index.html";
    my $index_file = 'index.html';
    my $index_entries = index_entry(\%genelist);
#(name,synonyme, 4 links, max-table wuerden reichen)
    my $index_vars = 
	{
	    genesofinterests => $index_entries
	};
    $template->process($index_file,$index_vars,$index_path) || die "Template process failed: ", $template->error(), "\n";
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
    (my $sample = $filetoparse) =~ s/\.csv//;
    my @samples = split(/_/,$sample);
    my %entries	    = %{$_[1]};
    print STDERR "Expression parsing $sample!\n";
    open (LIST,"<","$filetoparse");
    while(<LIST>){
	next if($_ =~ /^#/);
	my $line  = $_;
	my ($gene, $mb3, $mb7, $mb23, $eb3, $eb7, $eb23, $l3, $l7, $l23, $max, $mp3, $mp7, $mp23, $ep3, $ep7, $ep23) = split(/\t/,$line);
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
	push @{$entries{$gene}{$goto}{CUFFLINKS}{$sample}{$samples[1]}}, ($mb3, $mb7, $mb23);
	push @{$entries{$gene}{$goto}{CUFFLINKS}{$sample}{$samples[2]}}, ($eb3, $eb7, $eb23);
	push @{$entries{$gene}{$goto}{LOGEXPRESSION}{$sample}{LOG}}, ($l3, $l7, $l23);
	push @{$entries{$gene}{$goto}{MAX}{$sample}{PVAL}}, $max;
	push @{$entries{$gene}{$goto}{PEAKS}{$sample}{$samples[1]}}, ($mp3, $mp7, $mp23);
	push @{$entries{$gene}{$goto}{PEAKS}{$sample}{$samples[2]}}, ($ep3, $ep7, $ep23);
    }
    return (\%entries);
}

sub parse_comparison{
#HG19 RAE Mock
#HG19 RAE MARV
#HG19 RAE EBOV
##Gene basemean HuH7-3h basemean HuH7-7h basemean HuH7-23h basemean R06E-J-3h basemean R06E-J-7h basemean R06E-J-23h log2(fc) 3h log2(fc) 7h log2(fc) 23h padj	peak HuH7-3h peak HuH7-7h peak HuH7-23h	peak R06E-J-3h peak R06E-J-7h peak R06E-J-23h
    my $filetoparse = $_[0];
    (my $sample = $filetoparse) =~ s/\.csv//;
    my @samples = split(/_/,$sample);
    my %entries	    = %{$_[1]};
    print STDERR "Expression parsing $sample!\n";
    open (LIST,"<","$filetoparse");
    while(<LIST>){
	next if($_ =~ /^#/);
	my $line  = $_;
	my ($gene, $hg3, $hg7, $hg23, $rae3, $rae7, $rae23, $l3, $l7, $l23, $max, $hgp3, $hgp7, $hgp23, $raep3, $raep7, $raep23) = split(/\t/,$line);
	my $goto;
	if (defined $entries{$gene}{GOI}){
	    $goto = "GOI";
	}
	elsif(defined $entries{$gene}{APG}){
	    $goto = "APG";
	}
	else{
	    $goto = "COMPARISON";
	}
	push @{$entries{$gene}{$goto}{HB}{$sample}{$samples[0]."\_".$samples[2]}}, ($hg3, $hg7, $hg23);
	push @{$entries{$gene}{$goto}{HB}{$sample}{$samples[1]."\_".$samples[2]}}, ($rae3, $rae7, $rae23);
	push @{$entries{$gene}{$goto}{CLOGEXPRESSION}{$sample}{LOG}}, ($l3, $l7, $l23);
	push @{$entries{$gene}{$goto}{CMAX}{$sample}{PVAL}}, $max;
	push @{$entries{$gene}{$goto}{CPEAKS}{$sample}{$samples[0]."\_".$samples[2]}}, ($hgp3, $hgp7, $hgp23);
	push @{$entries{$gene}{$goto}{CPEAKS}{$sample}{$samples[1]."\_".$samples[2]}}, ($raep3, $raep7, $raep23);
    }
    return (\%entries);
}

sub parse_timepoints{
#rae_timepoints.csv
#hg19_timepoints.csv
##Gene	MOCK basemean 3h	MOCK basemean 7h	MOCK basemean 23h	EBOV basemean 3h	EBOV basemean 7h	EBOV basemean 23h	MARV basemean 3h	MARV basemean 7h	MARV basemean 23h	MOCK log2(fc) 3hvs7h	MOCK log2(fc) 3hvs23h	MOCK log2(fc) 7hvs23h	EBOV log2(fc) 3hvs7h	EBOV log2(fc) 3hvs23h	EBOV log2(fc) 7hvs23h	MARV log2(fc) 3hvs7h	MARV log2(fc) 3hvs23h	MARV log2(fc) 7hvs23h	padj(max(fc)	MOCK peak 3h	MOCK peak 7h	MOCK peak 23h	EBOV peak 3h	EBOV peak 7h	EBOV peak 23h	MARV peak 3h	MARV peak 7h	MARV peak 23h
    my $filetoparse = $_[0];
    (my $sample = $filetoparse) =~ s/\.csv//;
    my @samples = split(/_/,$sample);
    my %entries	    = %{$_[1]};
    print STDERR "Expression parsing $sample!\n";
    open (LIST,"<","$filetoparse");
    while(<LIST>){
	next if($_ =~ /^#/);
	my $line  = $_;
	my ($gene, $mb3, $mb7, $mb23, $eb3, $eb7, $eb23, $mv3, $mv7, $mv23, $ml3, $ml7, $ml23, $el3, $el7, $el23, $vl3, $vl7, $vl23, $max, $mp3, $mp7, $mp23, $ep3, $ep7, $ep23, $vp3, $vp7, $vp23) = split(/\t/,$line);
	my $goto;
	if (defined $entries{$gene}{GOI}){
	    $goto = "GOI";
	}
	elsif(defined $entries{$gene}{APG}){
	    $goto = "APG";
	}
	else{
	    $goto = "TIMEPOINTS";
	}
	push @{$entries{$gene}{$goto}{MEV}{$sample}{$samples[0]."\_".$samples[2]}}, ($mb3, $mb7, $mb23);
	push @{$entries{$gene}{$goto}{MEV}{$sample}{$samples[1]."\_".$samples[2]}}, ($eb3, $eb7, $eb23);
	push @{$entries{$gene}{$goto}{MEV}{$sample}{$samples[1]."\_".$samples[2]}}, ($mv3, $mv7, $mv23);
	push @{$entries{$gene}{$goto}{TLOGEXPRESSION}{$sample}{LOG}}, ($ml3, $ml7, $ml23);
	push @{$entries{$gene}{$goto}{TLOGEXPRESSION}{$sample}{LOG}}, ($el3, $el7, $el23);
	push @{$entries{$gene}{$goto}{TLOGEXPRESSION}{$sample}{LOG}}, ($vl3, $vl7, $vl23);
	push @{$entries{$gene}{$goto}{TMAX}{$sample}{PVAL}}, $max;
	push @{$entries{$gene}{$goto}{TPEAKS}{$sample}{$samples[0]."\_".$samples[2]}}, ($mp3, $mp7, $mp23);
	push @{$entries{$gene}{$goto}{TPEAKS}{$sample}{$samples[1]."\_".$samples[2]}}, ($ep3, $ep7, $ep23);
	push @{$entries{$gene}{$goto}{TPEAKS}{$sample}{$samples[1]."\_".$samples[2]}}, ($vp3, $vp7, $vp23);
    }
    return (\%entries);
}

sub read_tables{
    my $filetoparse = $_[0];
    my %entries	    = %{$_[1]};
    my @again;
    print STDERR "Parsing $filetoparse!\n";
    open (my $list,"<","$filetoparse");
    my @process = <$list>;

    if ($filetoparse =~ /goi/i){
	print STDERR "Processing GIO List!\n";
	while (@process) {
	    my $item = shift(@process);
#print STDERR $#process,"\n";
	    next unless ($item =~ /.goi./);
	    (my $line	  = $item) =~ s/,w+//g;
	    my @fields		  = split(/\;/,$line);
	    my $goi		  = $fields[0];
	    my $hacker		  = $fields[1];
	    next if ($hacker eq 'OPTIONAL');
	    my $gene		  = $fields[2];
	    next if (defined $entries{$gene}{GOI}{ID});

	    my $duplicate	  = $fields[3];   
#	    my @synonyms	  = split(/[,\s]+/,$fields[4]);
	    my @synonyms	  = (
		map {   
		    s/^\s+//;  # strip leading spaces
#                s/\s+$//;  # strip trailing spaces
		    $_         # return the modified string
		}
		split(/[,]+|[\s]{2,}|\t/,$fields[4])
		) if ($fields[4] ne '');
	    push @synonyms, $gene unless ($synonyms[0]);
	    if ($duplicate eq '' || $duplicate == 0 || ($duplicate && $fields[8] ne '')){
		my @pathways =	(
		    map {   
			s/^\s+//;  # strip leading spaces
#                s/\s+$//;  # strip trailing spaces
			$_         # return the modified string
		    }
		    split(/[,]+|[\s]{2,}|\t/,$fields[5])
		    ) if ($fields[5] ne '');
		push @pathways, 'Unknown' unless ($pathways[0]);
		my @literature =	(
		    map {   
			s/^\s+//;  # strip leading spaces
#                s/\s+$//;  # strip trailing spaces
			$_         # return the modified string
		    }
		    split(/[,]+|[\s]{2,}|\t/,$fields[6])
		    ) if ($fields[6] ne '');
		push @literature, 'Unknown' unless ($literature[0]);
		my $igvs	      = (split(/[:,\s\/]+/,$fields[7],2))[0] || '0';
		$igvs = 0 if ($igvs =~ /todo|NA/i);
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
		$ucscs = 0 if ($ucscs =~ /todo|NA/i);
		my $ucsc_conservation = $fields[18];
		my $hg_SNPs	      = $fields[19];
		my $rae_SNPs	      = $fields[20];
		my $sashimi	      = (split(/[:,\s\/]+/,$fields[21],2))[0] || '0';
		$sashimi = 0 if ($sashimi =~ /todo|NA/i);
		my $hg_intron	      = $fields[22];
		my $rae_intron	      = $fields[23];
		my $hg_5utr	      = $fields[24];
		my $hg_3utr	      = $fields[25];
		my $extra	      = $fields[26];
		my $notes	      = $fields[27];
		
		$entries{$gene}{GOI}{ID}  =	$goi;
		push @{$entries{$gene}{GOI}{SYNONYMS}}, @synonyms;
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
		    push @{$entries{$gene}{GOI}{UCSC}},"$goi\/snapshots/$goi\_ucsc.eps";
		}
		elsif ($ucscs == 0){
		    push @{$entries{$gene}{GOI}{UCSC}},"NONE";
		}
		else{
		    for (1..$ucscs){
			push @{$entries{$gene}{GOI}{UCSC}},"$goi\/snapshots/$goi\_ucsc$_\.eps";
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
		    $entries{$gene}{GOI} = $entries{$syn}{GOI} if (defined $entries{$syn}{GOI} && !defined $entries{$gene}{GOI});
		}
		push @process, $line unless (defined $entries{$gene}{GOI});
	    }
	}
    }
    elsif ($filetoparse =~ /apg/i){
	print STDERR "Processing APG List!\n";
	while (@process) {
	    my $item = shift(@process);
	    next unless ($item =~ /.apg./);
	    (my $line	  = $item) =~ s/,w+//g;
	    my @fields		  = split(/\;/,$line);
#Folder_ID;Hacker-Name;Gen;Duplicates?;Synonyms;Pathway;Literature;Screen IGV;Scale;Diff. Mock;Diff. Ebola;Diff. Marburg;Profile change;found?;Diff. Mock;Diff. Ebola;Diff. Marburg;Profile Change;Screen UCSC;UCSC Conservation;SNPs;SNPs;Sashimi plot;Intron transcripts;Intron transcripts;Upstream, 5'UTR;Downstream, 3'UTR;Extra Screens and why;Notes;Empty Entries?;Text Draft DONE;done.sh;Questions?;Done Lit Team / Done Questions;FRANZI FINAL CHECK;VERENA;STATE
#hg19.apg.00209;Stefanie W.;ABCA1;0;;;;4;0-376;ed;u17d;dd;0;1;00;00;00;low expression;1;11111;0;0;4;0;0;0;0;0;COOL:  significantly downregulated in EBOV23h (327-390-22);#VALUE!;1;1;0;;1;transport molecules across extra and intracellular membranes;CHECKED
	    my $apg		  = $fields[0];
	    my $hacker		  = $fields[1];
	    next if ($hacker eq 'OPTIONAL');
	    my $gene		  = $fields[2] =~ s/^\s+//g;;
	    next if (defined $entries{$gene}{APG}{ID});

	    my $duplicate	  = $fields[3];   
	    my @synonyms	  = (
		map {   
		    s/^\s+//;  # strip leading spaces
#                s/\s+$//;  # strip trailing spaces
		    $_         # return the modified string
		}
		split(/[,]+|[\s]{2,}|\t/,$fields[4])
		) if ($fields[4] ne '');
	    push @synonyms, $gene unless ($synonyms[0]);
	    if ($duplicate eq '' || $duplicate == 0){
		my @pathways =	(
		    map {   
			s/^\s+//;  # strip leading spaces
#                s/\s+$//;  # strip trailing spaces
			$_         # return the modified string
		    }
		    split(/[,]+|[\s]{2,}|\t/,$fields[5])
		    ) if ($fields[5] ne '');
		push @pathways, 'Unknown' unless ($pathways[0]);
		my @literature =	(
		    map {   
			s/^\s+//;  # strip leading spaces
#                s/\s+$//;  # strip trailing spaces
			$_         # return the modified string
		    }
		    split(/[,]+|[\s]{2,}|\t/,$fields[6])
		    ) if ($fields[6] ne '');
		push @literature, 'Unknown' unless ($literature[0]);
		my $igvs	      = (split(/[:,\s\/]+/,$fields[7],2))[0] || '0';
		$igvs = 0 if ($igvs =~ /todo|NA/i);
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
		my $ucscs	      = (split(/[:,\s\/]+/,$fields[18],2))[0] || '0';
		$ucscs = 0 if ($ucscs =~ /todo|NA/i);
		my $ucsc_conservation = $fields[19];
		my $hg_SNPs	      = $fields[20];
		my $rae_SNPs	      = $fields[21];
		my $sashimi	      = (split(/[:,\s\/]+/,$fields[22],2))[0] || '0';
		$sashimi = 0 if ($sashimi =~ /todo|NA/i);
		my $hg_intron	      = $fields[23];
		my $rae_intron	      = $fields[24];
		my $hg_5utr	      = $fields[25];
		my $hg_3utr	      = $fields[26];
		my $extra	      = $fields[27];
		my $notes	      = $fields[28];
		
		$entries{$gene}{APG}{ID}  = $apg;
		push @{$entries{$gene}{APG}{SYNONYMS}},  @synonyms;
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
		    push @{$entries{$gene}{APG}{UCSC}},"$apg\/snapshots/$apg\_ucsc.eps";
		}
		elsif ($ucscs == 0){
		    push @{$entries{$gene}{APG}{UCSC}},"NONE";
		}
		else{
		    for (1..$ucscs){
			push @{$entries{$gene}{APG}{UCSC}},"$apg\/snapshots/$apg\_ucsc$_\.eps";
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
		    $entries{$gene}{APG} = $entries{$syn}{APG} if (defined $entries{$syn}{APG} && !defined $entries{$gene}{APG});
		}
		push @process, $line unless (defined $entries{$gene}{APG});
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
    my $glist = shift;
    my %genelist = %{$glist};
    my $index_entry;
    foreach my $synonym (sort keys %genelist){
	my $goi = $genelist{$synonym};
	my $goilink = "$goi\.html";
	$index_entry .= "<tr><td><a href=\"$goilink\">$synonym</a></td></tr>";
    }
    return $index_entry;
}

sub image_entry{
    my $file = shift;
    my $dir = shift;
    my $odir = shift;
    my $image_entry;
    if($file eq 'NONE'){
        $image_entry = "none";
    }else{
        my @file = split ("/", $file);
        my $filename = $file[2];
        if($filename =~/.svg/){
            $filename =~ s/.svg/.png/;
            my $snapshotdir = join("/",$wdir,$dir,$file[0],$file[1]);
            my $imagelink = $wdir . "/" . $dir ."/". $file;
            my $thumblink = "./" ."thumbs/" . "$filename";
            `inkscape --file=$imagelink --export-width=150 --export-height=150 --without-gui --export-png=$thumblink` unless (-e $thumblink || !-e $imagelink);
            $image_entry = "<a href=\"$snapshotdir\"><img src=\"$thumblink\"></a>";
        }
        if($filename =~/.eps/){
            $filename =~ s/.eps/.png/;
            my $snapshotdir = join("/",$wdir,$dir,$file[0],$file[1]);
            my $imagelink = $wdir . "/" . $dir ."/". $file;
            my $thumblink = "./" ."thumbs/" . "$filename";
            `convert $imagelink -resize 150x150! $thumblink` unless (-e $thumblink || !-e $imagelink);
            $image_entry = "<a href=\"$snapshotdir\"><img src=\"$thumblink\"></a>";
        }
    }
    return $image_entry;
}

sub link_entry{
    my $file = shift;
    my $dir = shift;
    my $link_entry;
    if($file eq 'NONE'){
        $link_entry = "none";
    }else{
        my @file = split ("/", $file);
        my $filename = $file[2];
        my $snapshotdir = join("/",$wdir,$dir,$file[0],$file[1]);
        my $filelink = $wdir . "/" . $dir ."/". $file;
        $link_entry = "<a href=\"$filelink\">texfile</a>";
    }
    return $link_entry;
}
