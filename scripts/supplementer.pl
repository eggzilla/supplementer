#!/usr/bin/perl
### To use this Script replace all semicolons in fields with commas,
### then all line breaks in fields with nothing,
### then save as semicolon separated list and have fun parsing
### 
### Script supplementer.pl;
### Last changed Time-stamp: <2014-12-18 00:48:43 fall> by joerg

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
    elsif($file =~ /mock|ebov|marv/i && !($file =~ /hg19/i && $file =~ /rae/i)){
	%genes = %{parse_expression($file,\%genes)};
    }
    elsif($file =~ /timepoints/i){
	%genes = %{parse_timepoints($file,\%genes)};
    }
    elsif($file =~ /hg19/i && $file =~ /rae/i){
	%genes = %{parse_comparison($file,\%genes)};
    }
    elsif($file =~ /deseq/i){
	%genes = %{parse_deseq($file,\%genes)};
    }
    elsif($file =~ /ymax.norm/i){
	%genes = %{parse_ymax($file,\%genes)};
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
    my @parseit = ('GOI', 'APG', 'EXPRESSION', 'TIMEPOINTS' , 'COMPARISON', 'DESEQ');
#    my @sorted_genes = 	map { $_->[0] } sort { $a->[1] cmp $b->[1] } map { [ $_, uc($_) ] } keys %gois;
    my $index_entries;
    my %peaks;
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
	    my (@samp, @condi, @deg, @max, @maxl, @fold) = ();
	    foreach my $sample (sort {lc($a) cmp lc($b)} keys %{$gois{$gene}{$from}{CUFFLINKS}} ){
		my $fold_change = 'NA';
		$fold_change = join(" | ",@{$gois{$gene}{$from}{LOGEXPRESSION}{$sample}{LOG}}) if (defined $gois{$gene}{$from}{LOGEXPRESSION}{$sample}{LOG});
		push @fold, $fold_change;
		push @samp, $sample;
		foreach my $condition (sort {lc($a) cmp lc($b)} keys %{$gois{$gene}{$from}{CUFFLINKS}{$sample}} ){
		    my $cufflinks = 'NA';
		    $cufflinks = join(" | ",@{$gois{$gene}{$from}{CUFFLINKS}{$sample}{$condition}}) if (defined $gois{$gene}{$from}{CUFFLINKS}{$sample}{$condition});
		    @{$gois{$gene}{$from}{PEAKS}{$sample}{$condition}} = grep /\S/, @{$gois{$gene}{$from}{PEAKS}{$sample}{$condition}} if ($gois{$gene}{$from}{PEAKS}{$sample}{$condition}); ## get rid of empty entries
		    my $maxy = 'NA';
		    $maxy = join(" | ",@{$gois{$gene}{$from}{PEAKS}{$sample}{$condition}}) if ($gois{$gene}{$from}{PEAKS}{$sample}{$condition});
		    @{$peaks{$gene}{$from}{$condition}} = @{$gois{$gene}{$from}{PEAKS}{$sample}{$condition}} if ($gois{$gene}{$from}{PEAKS}{$sample}{$condition} && $sample =~ /hg19/);
		    push @deg, $cufflinks;
		    push @max, $maxy;
		    push @condi, $condition;
		}
	    }
### Parse peaks separate for each condition, independent of sample to get rid of reocurring values
	    foreach my $condition ( sort {lc($a) cmp lc($b)} keys %{$peaks{$gene}{$from}} ){
		my $maxy = 'NA';
		$maxy = join(",",@{$peaks{$gene}{$from}{$condition}}) if ($peaks{$gene}{$from}{$condition});
		push @maxl, $maxy;
	    }
	    my $peak;
	    my @ranges = (2,0,1); ## Resort order for index.html page
	    for (@ranges){
		$peak .= $maxl[$_]."," if ($maxl[$_]);
	    }
	    $peak = 'NA' unless ($peak && $peak !~ /NA/i);

### Parse Comparison
	    my (@csamp, @ccondi, @cdeg, @cmax, @cfold) = ();
	    foreach my $sample (sort {lc($a) cmp lc($b)} keys %{$gois{$gene}{$from}{HB}} ){
		my $cfold_change = 'NA';
		$cfold_change = join(" | ",@{$gois{$gene}{$from}{CLOGEXPRESSION}{$sample}{LOG}}) if(defined $gois{$gene}{$from}{CLOGEXPRESSION}{$sample}{LOG});
		push @cfold, $cfold_change;
		push @csamp, $sample;
		foreach my $condition (sort {lc($a) cmp lc($b)} keys %{$gois{$gene}{$from}{HB}{$sample}} ){
		    my $cufflinks = 'NA';
		    $cufflinks = join(" | ",@{$gois{$gene}{$from}{HB}{$sample}{$condition}}) if (defined $gois{$gene}{$from}{HB}{$sample}{$condition});
		    push @cdeg, $cufflinks;	    
		    push @ccondi, $condition;
		}
		foreach my $condition (sort {lc($a) cmp lc($b)} keys %{$gois{$gene}{$from}{CPEAKS}{$sample}} ){
		    @{$gois{$gene}{$from}{CPEAKS}{$sample}{$condition}} = grep /\S/,@{$gois{$gene}{$from}{CPEAKS}{$sample}{$condition}} if ($gois{$gene}{$from}{CPEAKS}{$sample}{$condition});
		    my $maxy = 'NA';
		    $maxy = join(" | ",@{$gois{$gene}{$from}{CPEAKS}{$sample}{$condition}}) if (defined $gois{$gene}{$from}{CPEAKS}{$sample}{$condition});
		    push @cmax, $maxy;
		}
	    }

### Parse Timepoints
	    my (@tsamp, @tcondi, @tdeg, @tmax, @tfold) = ();
	    foreach my $sample (keys %{$gois{$gene}{$from}{MEV}} ){
		push @tsamp, $sample;
		foreach my $condition (sort {lc($a) cmp lc($b)} keys %{$gois{$gene}{$from}{MEV}{$sample}} ){
		    my $cufflinks = 'NA';
		    $cufflinks = join(" | ",@{$gois{$gene}{$from}{MEV}{$sample}{$condition}}) if (defined $gois{$gene}{$from}{MEV}{$sample}{$condition});
		    @{$gois{$gene}{$from}{TPEAKS}{$sample}{$condition}} = grep /\S/, @{$gois{$gene}{$from}{TPEAKS}{$sample}{$condition}} if ($gois{$gene}{$from}{TPEAKS}{$sample}{$condition});
		    my $maxy = 'NA';
		    $maxy = join(" | ",@{$gois{$gene}{$from}{TPEAKS}{$sample}{$condition}}) if (defined $gois{$gene}{$from}{TPEAKS}{$sample}{$condition});
		    my $fold_change = 'NA';
		    $fold_change = join(" | ",@{$gois{$gene}{$from}{TLOGEXPRESSION}{$sample}{$condition}}) if (defined $gois{$gene}{$from}{TLOGEXPRESSION}{$sample}{$condition});
		    push @tdeg, $cufflinks;
		    push @tmax, $maxy;
		    push @tfold, $fold_change;
		    push @tcondi, $condition;
		}
	    }

### Parse DESeq
	    my (@dsamp, @dcondi, @ddeg, @dfold, @defold) = ();
	    foreach my $sample (sort {lc($a) cmp lc($b)} keys %{$gois{$gene}{$from}{DLOGEXPRESSION}} ){
		push @dsamp, $sample;
		foreach my $condition (sort {lc($a) cmp lc($b)} keys %{$gois{$gene}{$from}{DLOGEXPRESSION}{$sample}} ){
		    my $cufflinks = 'NA';
		    $cufflinks = join(" | ",@{$gois{$gene}{$from}{DE}{$sample}{$condition}}) if (defined $gois{$gene}{$from}{DE}{$sample}{$condition});
		    my $fold_change = 'NA';
		    $fold_change = join(" | ",@{$gois{$gene}{$from}{DLOGEXPRESSION}{$sample}{$condition}}) if (defined $gois{$gene}{$from}{DLOGEXPRESSION}{$sample}{$condition});
		    push @ddeg, $cufflinks;
		    push @dfold, $fold_change;
		    push @dcondi, $condition;
		}
		foreach my $condition (sort {lc($a) cmp lc($b)} keys %{$gois{$gene}{$from}{DELOGEXPRESSION}{$sample}} ){
		    my $defold_change = 'NA';
		    $defold_change = join(" | ",@{$gois{$gene}{$from}{DELOGEXPRESSION}{$sample}{$condition}}) if (defined $gois{$gene}{$from}{DELOGEXPRESSION}{$sample}{$condition});
		    push @defold, $defold_change;
		}
	    }
	    @ddeg = ($deg[0],$deg[1],$deg[3]) unless (($ddeg[0] && $ddeg[0] ne '') || !$deg[0]);
	    
            my $igv = 'NONE';
	    $igv = image_entry(${$gois{$gene}{$from}{IGV}}[0],$dir,$odir) if ($from eq 'GOI' || $from eq 'APG');
            my $sashimi =  'NONE';
	    $sashimi = image_entry(${$gois{$gene}{$from}{SASHIMI}}[0],$dir,$odir) if ($from eq 'GOI' || $from eq 'APG');
            my $ucsc =  'NONE';
	    $ucsc = image_entry(${$gois{$gene}{$from}{UCSC}}[0],$dir,$odir) if ($from eq 'GOI' || $from eq 'APG');
            my $tex_link =  'NONE';
	    $tex_link = link_entry($gois{$gene}{$from}{TEX},$dir) if ($from eq 'GOI' || $from eq 'APG');
	    my $syn = 'UNKNOWN';
	    ($syn = join(",",@{$gois{$gene}{$from}{SYNONYMS}})) =~ s/, /,/g if ($from eq 'GOI' || $from eq 'APG');
	    my $goi_link = goi_link($gene,$gois{$gene}{$from}{ID});
	    $index_entries .= index_entry_detailed($template_path,$goi_link,$syn,$gois{$gene}{$from}{ID},$tex_link,$igv,$sashimi,$ucsc,$peak);
            my $texcontent = 'NA';
	    $texcontent = tex_content($wdir,$dir,$gois{$gene}{$from}{TEX}) if ($from eq 'GOI' || $from eq 'APG');
	    my $add = 'NA';
	    $add = additional_plot_entry($gois{$gene}{$from}{EXTRA},${$gois{$gene}{$from}{IGV}}[0],$dir) if ($from eq 'GOI' || $from eq 'APG');

	    my $goi_vars = 
	    {   
		name		  => $gois{$gene}{$from}{NAME} || 'NA',
		synonyms	  => $syn || 'NA',
		goiid		  => $gois{$gene}{$from}{ID} || 'NA',
		maxy		  => $peak || 'NA',
		textxt		  => $tex_link || 'NA',
		igv		  => $igv || 'NA',
		sashimi		  => $sashimi || 'NA',
		ucsc		  => $ucsc || 'NA',
		additionalplots   => $add || 'NA',
                texcontent        => $texcontent || 'NA',
##sample1
		sample		  => $samp[0] || 'NA',
##sample 1 condition1
		condone		  => $condi[0] || 'NA',
		maxy		  => $max[0] || 'NA',
		cufflinks	  => $deg[0] || 'NA',
		fold		  => $fold[0] || 'NA',
##sample1 condition2
		condtwo		  => $condi[1] || 'NA',
		maxy_two	  => $max[1] || 'NA',
		cufflinks_two     => $deg[1] || 'NA',
##sample2
		sampleone	  => $samp[1] || 'NA',
##sample 2 condition1
		condoone	  => $condi[2] || 'NA',
		maxyo		  => $max[2] || 'NA',
		cufflinkso	  => $deg[2] || 'NA',
		foldo		  => $fold[1] || 'NA',
##sample 2 condition2
		condotwo	  => $condi[3] || 'NA',
		maxyo_two	  => $max[3] || 'NA',
		cufflinkso_two    => $deg[3] || 'NA',
##sample 3
		sampletwo	  => $samp[2] || 'NA',
##sample 3 condition1
		condtone	  => $condi[4] || 'NA',
		maxyt		  => $max[4] || 'NA',
		cufflinkst	  => $deg[4] || 'NA',
		foldt		  => $fold[2] || 'NA',    
##sample 3 condition2
		condttwo	  => $condi[5] || 'NA',
		maxyt_two	  => $max[5] || 'NA',
		cufflinkst_two    => $deg[5] || 'NA',
###Comparison
##sample1
		csample		  => $csamp[0] || 'NA',
##sample 1 condition1
		ccondone	  => $ccondi[0] || 'NA',
		cmaxy		  => $cmax[0] || 'NA',
		ccufflinks	  => $cdeg[0] || 'NA',
		cfold             => $cfold[0] || 'NA',
##sample1 condition2
		ccondtwo	  => $ccondi[1] || 'NA',
		cmaxy_two	  => $cmax[1] || 'NA',
		ccufflinks_two    => $cdeg[1] || 'NA',
##sample2
		csampleone	  => $csamp[1] || 'NA',
##sample 2 condition1
		ccondoone	  => $ccondi[2] || 'NA',
		cmaxyo		  => $cmax[2] || 'NA',
		ccufflinkso	  => $cdeg[2] || 'NA',
		cfoldo            => $cfold[1] || 'NA',
##sample 2 condition2
		ccondotwo	  => $ccondi[3] || 'NA',
		cmaxyo_two	  => $cmax[3] || 'NA',
		ccufflinkso_two   => $cdeg[3] || 'NA',
### Timepoints
##sample1
		tsample		  => $tsamp[0] || 'NA',
##sample 1 condition1
		tcondone	  => $tcondi[0] || 'NA',
		tmaxy		  => $tmax[0] || 'NA',
		tcufflinks	  => $tdeg[0] || 'NA',
		tfold             => $tfold[0] || 'NA',
##sample1 condition2
		tcondtwo	  => $tcondi[1] || 'NA',
		tmaxy_two	  => $tmax[1] || 'NA',
		tcufflinks_two    => $tdeg[1] || 'NA',
		tfold_two         => $tfold[1] || 'NA',
##sample1 condition3
		tcondthree	  => $tcondi[2] || 'NA',
		tmaxy_three	  => $tmax[2] || 'NA',
		tcufflinks_three  => $tdeg[2] || 'NA',
		tfold_three       => $tfold[2] || 'NA',
### DESeq
##sample1
		dsample		  => $dsamp[0] || 'NA',
##sample 1 condition1
		dcondone	  => $dcondi[0] || 'NA',
#		dmaxy		  => $dmax[0],
		dcufflinks	  => $ddeg[0] || 'NA',
		dfold             => $dfold[0] || 'NA',
		defold            => $defold[0] || 'NA',
##sample1 condition2
		dcondtwo	  => $dcondi[1] || 'NA',
#		dmaxy_two	  => $dmax[1],
		dcufflinks_two    => $ddeg[1] || 'NA',
		dfold_two         => $dfold[1] || 'NA',
		defold_two        => $defold[1] || 'NA',
##sample1 condition3
		dcondthree	  => $dcondi[2] || 'NA',
#		dmaxy_three	  => $dmax[2],
		dcufflinks_three  => $ddeg[2] || 'NA',
		dfold_three       => $dfold[2] || 'NA',
		defold_three      => $defold[2] || 'NA'
	    };

	    $template->process($goi_file,$goi_vars,$goi_path) || die "Template process failed: ", $template->error(), "\n";	
	}
    }
#construct index.hmtl
    my $index_path = $html_destination_path. "/index.html";
    my $index_file = 'index.html';
    my $index_vars = 
    {
	genesofinterests => $index_entries
    };
    $template->process($index_file,$index_vars,$index_path) || die "Template process failed: ", $template->error(), "\n";
    chdir($wdir) or die "$!";
}

sub tex_content{
    my $wdir = shift;
    my $dir = shift;
    my $file = shift;
    my $filetoparse = $wdir . "/" . $dir ."/". $file;
    my @description;
    my $read=0;
    if (-e $filetoparse){
	open (LIST,"<","$filetoparse") || die "$!";
	while(<LIST>){
	    chomp(my $line  = $_);
	    next if ($line =~ /genetoclevel|put your text|^$/);
	    push @description, $line unless ($line=~/^%%%%%%%%%%%%%%%%%%%%%/);	
	    last if ($line=~/^%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/);
	}
    }
    else{
	push @description, 'No TeX file found!';
    }
    my $descript = join("\n",@description);
    return $descript;
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
    print STDERR "Expression parsing $sample @samples!\n";
    open (LIST,"<","$filetoparse") || die "$!\n";
    while(<LIST>){
	next if($_ =~ /^#/);
	my $line  = $_;
	my ($gene, $mb3, $mb7, $mb23, $eb3, $eb7, $eb23, $l3, $l7, $l23, $max, $mp3, $mp7, $mp23, $ep3, $ep7, $ep23) = split(/\t/,$line);
	my $goto;
	if (defined $entries{$gene}{GOI}{ID}){
	    $goto = "GOI";
	}
	elsif(defined $entries{$gene}{APG}{ID}){
	    $goto = "APG";
	}
	elsif((keys %{$entries{$gene}})[0]){   
	    $goto = (keys %{$entries{$gene}})[0];
	}
	else{
	    $goto = "EXPRESSION";
	}
	$entries{$gene}{$goto}{ID} = $gene if ($goto eq 'EXPRESSION' && $gene ne '');
	$entries{$gene}{$goto}{NAME} = $gene if ($goto eq 'EXPRESSION' && $gene ne '');
	push @{$entries{$gene}{$goto}{CUFFLINKS}{$sample}{$samples[1]}}, ($mb3, $mb7, $mb23);
	push @{$entries{$gene}{$goto}{CUFFLINKS}{$sample}{$samples[2]}}, ($eb3, $eb7, $eb23);
	push @{$entries{$gene}{$goto}{LOGEXPRESSION}{$sample}{LOG}}, ($l3, $l7, $l23);
	push @{$entries{$gene}{$goto}{MAX}{$sample}{PVAL}}, $max;
	push @{$entries{$gene}{$goto}{PEAKS}{$sample}{$samples[1]}}, ($mp3, $mp7, $mp23);
	push @{$entries{$gene}{$goto}{PEAKS}{$sample}{$samples[2]}}, ($ep3, $ep7, $ep23);

	my @syno;
	push @syno, @{$entries{$gene}{$goto}{SYNONYMS}} if (defined $entries{$gene}{$goto}{SYNONYMS});
	foreach my $syn (@syno){
	    next if ($syn eq $gene);
	    $entries{$syn}{$goto}{LOGEXPRESSION}=$entries{$gene}{$goto}{LOGEXPRESSION} if (defined $entries{$gene}{$goto}{LOGEXPRESSION} && !defined $entries{$syn}{$goto}{LOGEXPRESSION});
	    $entries{$syn}{$goto}{CUFFLINKS}=$entries{$gene}{$goto}{CUFFLINKS} if (defined $entries{$gene}{$goto}{CUFFLINKS} && !defined $entries{$syn}{$goto}{CUFFLINKS});
	    $entries{$syn}{$goto}{PEAKS}=$entries{$gene}{$goto}{PEAKS} if (defined $entries{$gene}{$goto}{PEAKS} && !defined $entries{$syn}{$goto}{PEAKS});
	} 
    }
    return (\%entries);
}

sub parse_ymax{
#hg19.ymax.norm
    my $filetoparse = $_[0];
    my %entries	    = %{$_[1]};
    print STDERR "Ymax parsing $filetoparse!\n";
    open (LIST,"<","$filetoparse") || die "$!\n";
    while(<LIST>){
	next if($_ =~ /^#/);
	my $line  = $_;
	my ($gene, $mp3, $mp7, $mp23, $ep3, $ep7, $ep23, $vp3, $vp7, $vp23) = split(/\t/,$line);
	my $goto;
	if (defined $entries{$gene}{GOI}{ID}){
	    $goto = "GOI";
	}
	elsif(defined $entries{$gene}{APG}{ID}){
	    $goto = "APG";
	}
	else{
	    next;
	}
	my @samp = ('hg19_mock_ebov','hg19_mock_marv','hg19_ebov_marv');
	my @orgs = ('mock','ebov','marv');
	foreach my $sample (@samp){
	    foreach my $org (@orgs){
		push @{$entries{$gene}{$goto}{PEAKS}{$sample}{$org}}, ($mp3, $mp7, $mp23) if ($org eq 'mock');
		push @{$entries{$gene}{$goto}{PEAKS}{$sample}{$org}}, ($ep3, $ep7, $ep23) if ($org eq 'ebov');
		push @{$entries{$gene}{$goto}{PEAKS}{$sample}{$org}}, ($vp3, $vp7, $vp23) if ($org eq 'marv');
	    }
	}
	my @syno;
	push @syno, @{$entries{$gene}{$goto}{SYNONYMS}} if (defined $entries{$gene}{$goto}{SYNONYMS});
	foreach my $syn (@syno){
	    next if ($syn eq $gene);
	    $entries{$syn}{$goto}{PEAKS}=$entries{$gene}{$goto}{PEAKS} if (defined $entries{$gene}{$goto}{PEAKS} && !defined $entries{$syn}{$goto}{PEAKS});
	}

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
    print STDERR "Comparison parsing $sample!\n";
    open (LIST,"<","$filetoparse") || die "$!\n";
    while(<LIST>){
	next if($_ =~ /^#/);
	my $line  = $_;
	my ($gene, $hg3, $hg7, $hg23, $rae3, $rae7, $rae23, $l3, $l7, $l23, $max, $hgp3, $hgp7, $hgp23, $raep3, $raep7, $raep23) = split(/\t/,$line);
	my $goto;
	if (defined $entries{$gene}{GOI}{ID}){
	    $goto = "GOI";
	}
	elsif(defined $entries{$gene}{APG}{ID}){
	    $goto = "APG";
	}
	elsif((keys %{$entries{$gene}})[0]){
	    $goto = (keys %{$entries{$gene}})[0];
	}
	else{
	    $goto = "COMPARISON";
	}
	$entries{$gene}{$goto}{ID} = $gene if ($goto eq 'COMPARISON' && $gene ne '');
	$entries{$gene}{$goto}{NAME} = $gene if ($goto eq 'COMPARISON' && $gene ne '');
	push @{$entries{$gene}{$goto}{HB}{$sample}{$samples[0]."\_".$samples[2]}}, ($hg3, $hg7, $hg23);
	push @{$entries{$gene}{$goto}{HB}{$sample}{$samples[1]."\_".$samples[2]}}, ($rae3, $rae7, $rae23);
	push @{$entries{$gene}{$goto}{CLOGEXPRESSION}{$sample}{LOG}}, ($l3, $l7, $l23);
	push @{$entries{$gene}{$goto}{CMAX}{$sample}{PVAL}}, $max;
	push @{$entries{$gene}{$goto}{CPEAKS}{$sample}{$samples[0]."\_".$samples[2]}}, ($hgp3, $hgp7, $hgp23);
	push @{$entries{$gene}{$goto}{CPEAKS}{$sample}{$samples[1]."\_".$samples[2]}}, ($raep3, $raep7, $raep23);
	
	my @syno;
	push @syno, @{$entries{$gene}{$goto}{SYNONYMS}} if (defined $entries{$gene}{$goto}{SYNONYMS});
	foreach my $syn (@syno){
	    next if ($syn eq $gene);
	    $entries{$syn}{$goto}=$entries{$gene}{$goto};		    
	} 
    }
    return (\%entries);
}

sub parse_timepoints{
#rae_timepoints.csv
#hg19_timepoints.csv
##Gene	MOCK basemean 3h	MOCK basemean 7h	MOCK basemean 23h	EBOV basemean 3h	EBOV basemean 7h	EBOV basemean 23h	MARV basemean 3h	MARV basemean 7h	MARV basemean 23h	MOCK log2(fc) 3hvs7h	MOCK log2(fc) 3hvs23h	MOCK log2(fc) 7hvs23h	EBOV log2(fc) 3hvs7h	EBOV log2(fc) 3hvs23h	EBOV log2(fc) 7hvs23h	MARV log2(fc) 3hvs7h	MARV log2(fc) 3hvs23h	MARV log2(fc) 7hvs23h	padj(max(fc)	MOCK peak 3h	MOCK peak 7h	MOCK peak 23h	EBOV peak 3h	EBOV peak 7h	EBOV peak 23h	MARV peak 3h	MARV peak 7h	MARV peak 23h
    my $filetoparse = $_[0];
    (my $sample = $filetoparse) =~ s/\.csv//;
    my @samples = ('mock','ebov','marv');
    my %entries	    = %{$_[1]};
    print STDERR "Timepoint parsing $sample!\n";
    open (LIST,"<","$filetoparse") || die "$!\n";
    while(<LIST>){
	next if($_ =~ /^#/);
	my $line  = $_;
	my ($gene, $mb3, $mb7, $mb23, $eb3, $eb7, $eb23, $mv3, $mv7, $mv23, $ml3, $ml7, $ml23, $el3, $el7, $el23, $vl3, $vl7, $vl23, $max, $mp3, $mp7, $mp23, $ep3, $ep7, $ep23, $vp3, $vp7, $vp23) = split(/\t/,$line);
	my $goto;
	if (defined $entries{$gene}{GOI}{ID}){
	    $goto = "GOI";
	}
	elsif(defined $entries{$gene}{APG}{ID}){
	    $goto = "APG";
	}
	elsif((keys %{$entries{$gene}})[0]){
	    $goto = (keys %{$entries{$gene}})[0];
	}
	else{
	    $goto = "TIMEPOINTS";
	}
	$entries{$gene}{$goto}{ID} = $gene if ($goto eq 'TIMEPOINTS' && $gene ne '');
	$entries{$gene}{$goto}{NAME} = $gene if ($goto eq 'TIMEPOINTS' && $gene ne '');
	push @{$entries{$gene}{$goto}{MEV}{$sample}{$samples[0]}}, ($mb3, $mb7, $mb23);
	push @{$entries{$gene}{$goto}{MEV}{$sample}{$samples[1]}}, ($eb3, $eb7, $eb23);
	push @{$entries{$gene}{$goto}{MEV}{$sample}{$samples[2]}}, ($mv3, $mv7, $mv23);
	push @{$entries{$gene}{$goto}{TLOGEXPRESSION}{$sample}{$samples[0]}}, ($ml3, $ml7, $ml23);
	push @{$entries{$gene}{$goto}{TLOGEXPRESSION}{$sample}{$samples[1]}}, ($el3, $el7, $el23);
	push @{$entries{$gene}{$goto}{TLOGEXPRESSION}{$sample}{$samples[2]}}, ($vl3, $vl7, $vl23);
	push @{$entries{$gene}{$goto}{TMAX}{$sample}{PVAL}}, $max;
	push @{$entries{$gene}{$goto}{TPEAKS}{$sample}{$samples[0]}}, ($mp3, $mp7, $mp23);
	push @{$entries{$gene}{$goto}{TPEAKS}{$sample}{$samples[1]}}, ($ep3, $ep7, $ep23);
	push @{$entries{$gene}{$goto}{TPEAKS}{$sample}{$samples[2]}}, ($vp3, $vp7, $vp23);

	my @syno;
	push @syno, @{$entries{$gene}{$goto}{SYNONYMS}} if (defined $entries{$gene}{$goto}{SYNONYMS});
	foreach my $syn (@syno){
	    next if ($syn eq $gene);
	    $entries{$syn}{$goto}=$entries{$gene}{$goto};		    
	} 
    }
    return (\%entries);
}

sub parse_deseq{
##Gene   Mock_3h Mock_7h Mock_23h        EBOV_3h EBOV_7h EBOV_23h        MARV_3h MARV_7h MARV_23h        FC(Mock-3h_vs_Mock-7h)  FC(Mock-3h_vs_Mock-23h) FC(Mock-7h_vs_Mock-23h) FC(EBOV-3h_vs_EBOV-7h)  FC(EBOV-3h_vs_EBOV-23h) FC(EBOV-7h_vs_EBOV-23h) FC(MARV-3h_vs_MARV-7h)  FC(MARV-3h_vs_MARV-23h) FC(MARV-7h_vs_MARV-23h) 
#FC(Mock-3h_vs_EBOV-3h)  FC(Mock-3h_vs_MARV-3h)  FC(EBOV-3h_vs_MARV-3h)  FC(Mock-7h_vs_EBOV-7h)  FC(Mock-7h_vs_MARV-7h)  FC(EBOV-7h_vs_MARV-7h)  FC(Mock-23h_vs_EBOV-23h)        FC(Mock-23h_vs_MARV-23h)        FC(EBOV-23h_vs_MARV-23h)        PADJ(Mock-3h_vs_Mock-7h)        PADJ(Mock-3h_vs_Mock-23h)       PADJ(Mock-7h_vs_Mock-23h)       PADJ(EBOV-3h_vs_EBOV-7h)        PADJ(EBOV-3h_vs_EBOV-23h)       PADJ(EBOV-7h_vs_EBOV-23h)       PADJ(MARV-3h_vs_MARV-7h)        PADJ(MARV-3h_vs_MARV-23h)       PADJ(MARV-7h_vs_MARV-23h)       PADJ(Mock-3h_vs_EBOV-3h)        PADJ(Mock-3h_vs_MARV-3h)        PADJ(EBOV-3h_vs_MARV-3h)        PADJ(Mock-7h_vs_EBOV-7h)        PADJ(Mock-7h_vs_MARV-7h)        PADJ(EBOV-7h_vs_MARV-7h)        PADJ(Mock-23h_vs_EBOV-23h)      PADJ(Mock-23h_vs_MARV-23h)      PADJ(EBOV-23h_vs_MARV-23h)
    my $filetoparse = $_[0];
    (my $sample = $filetoparse) =~ s/\.csv//;
    my @samples = split(/\_/,$sample);
    my %entries	    = %{$_[1]};
    print STDERR "Timepoint parsing $sample!\n";
    open (LIST,"<","$filetoparse") || die "$!\n";
    while(<LIST>){
	next if($_ =~ /^#/);
	my $line  = $_;
	my ( $gene, $mb3, $mb7, $mb23, $eb3, $eb7, $eb23, $vb3, $vb7, $vb23, 
	     $ml37, $ml323, $ml723, $el37, $el323, $el723, $vl37, $vl323, $vl723,
	     $me3, $mv3, $ev3, $me7, $mv7, $ev7, $me23, $mv23, $ev23 ) = split(/\t/,$line);
#	next if ($mb3 =~ /NaN|Infinity/ || $el37 =~ /NaN|Infinity/);
	my $goto;
	if (defined $entries{$gene}{GOI}){
	    $goto = "GOI";
	}
	elsif(defined $entries{$gene}{APG}){
	    $goto = "APG";
	}
	elsif((keys %{$entries{$gene}})[0]){
	    $goto = (keys %{$entries{$gene}})[0];
	}
	else{
	    $goto = "DESEQ";
	}
	$entries{$gene}{$goto}{ID} = $gene if ($goto eq 'DESEQ' && $gene ne '');
	$entries{$gene}{$goto}{NAME} = $gene if ($goto eq 'DESEQ' && $gene ne '');
	my $sampled;
	push @{$entries{$gene}{$goto}{DE}{$sample}{mock}}, ($mb3, $mb7, $mb23);
	push @{$entries{$gene}{$goto}{DE}{$sample}{ebov}}, ($eb3, $eb7, $eb23);
	push @{$entries{$gene}{$goto}{DE}{$sample}{marv}}, ($vb3, $vb7, $vb23);
	if (!defined $entries{$gene}{$goto}{CUFFLINKS}{hg19_mock_ebov}{mock}){
	    $sampled = 'hg19_mock_ebov';
	    push @{$entries{$gene}{$goto}{CUFFLINKS}{$sampled}{mock}}, ($mb3, $mb7, $mb23);
	    push @{$entries{$gene}{$goto}{CUFFLINKS}{$sampled}{ebov}}, ($eb3, $eb7, $eb23);
	}
	if (!defined $entries{$gene}{$goto}{CUFFLINKS}{hg19_mock_marv}{mock}){
	    $sampled = 'hg19_mock_marv';
	    push @{$entries{$gene}{$goto}{CUFFLINKS}{$sampled}{mock}}, ($mb3, $mb7, $mb23);
	    push @{$entries{$gene}{$goto}{CUFFLINKS}{$sampled}{marv}}, ($vb3, $vb7, $vb23);
	}
	if (!defined $entries{$gene}{$goto}{CUFFLINKS}{hg19_ebov_marv}{ebov}){
	    $sampled = 'hg19_ebov_marv';
	    push @{$entries{$gene}{$goto}{CUFFLINKS}{$sampled}{ebov}}, ($eb3, $eb7, $eb23);
	    push @{$entries{$gene}{$goto}{CUFFLINKS}{$sampled}{marv}}, ($vb3, $vb7, $vb23);
	}
	push @{$entries{$gene}{$goto}{DLOGEXPRESSION}{$sample}{mock}}, ($ml37, $ml323, $ml723);
	push @{$entries{$gene}{$goto}{DLOGEXPRESSION}{$sample}{ebov}}, ($el37, $el323, $el723);
	push @{$entries{$gene}{$goto}{DLOGEXPRESSION}{$sample}{marv}}, ($vl37, $vl323, $vl723);
	push @{$entries{$gene}{$goto}{DELOGEXPRESSION}{$sample}{DE3}}, ($me3, $mv3, $ev3);
	push @{$entries{$gene}{$goto}{DELOGEXPRESSION}{$sample}{DE7}}, ($me7, $mv7, $ev7);
	push @{$entries{$gene}{$goto}{DELOGEXPRESSION}{$sample}{DE23}}, ($me23, $mv23, $ev23);
	
	my @syno;
	push @syno, @{$entries{$gene}{$goto}{SYNONYMS}} if (defined $entries{$gene}{$goto}{SYNONYMS});
	foreach my $syn (@syno){
	    next if ($syn eq $gene);
	    $entries{$syn}{$goto}=$entries{$gene}{$goto};		    
	} 
    }
    return (\%entries);
}

sub read_tables{
    my $filetoparse = $_[0];
    my %entries	    = %{$_[1]};
    my @again;
    print STDERR "Parsing $filetoparse!\n";
    open (my $list,"<","$filetoparse") || die "$!\n";
    my @process = <$list>;

    if ($filetoparse =~ /goi/i){
	print STDERR "Processing GIO List!\n";
	while (@process) {
	    my $item = shift(@process);
	    next unless ($item =~ /.goi./);
	    my $to = 'GOI';
	    (my $line	  = $item) =~ s/,w+//g;
	    my @fields		  = split(/\;/,$line);
	    my $goi		  = $fields[0];
	    my $hacker		  = $fields[1];
	    next if ($hacker eq 'OPTIONAL' || $hacker eq '');
	    (my $gene		  = $fields[2]) =~ s/^\s+//g;
	    $gene = $goi unless $gene;;
	    next if (defined $entries{$gene}{$to}{ID});
#	    $to = 'APG' if (exists $entries{$gene}{APG});
	    my $duplicate = $fields[3];   
	    my @synonyms  = ();
	    @synonyms	 = (
		map {   
		    s/^\s+//;  # strip leading spaces
		    s/\s+$//;  # strip trailing spaces
		    $_         # return the modified string
		}
		split(/[,]+|[\s]{2,}|\t/,$fields[4])
		) if ($fields[4] ne '');
	    push @synonyms, $gene unless ($gene =~ /.goi./ || $gene eq '');
	    @synonyms = grep /\S/, @synonyms; ## get rid of empty entries
	    if ($duplicate eq '' || $duplicate == 0 || ($duplicate && $fields[8] ne '')){
		my @pathways = ();
		@pathways    = (
		    map {   
			s/^\s+//;  # strip leading spaces
			s/\s+$//;  # strip trailing spaces
			$_         # return the modified string
		    }
		    split(/[,]+|[\s]{2,}|\t/,$fields[5])
		    ) if ($fields[5] ne '');
		push @pathways, 'Unknown' unless ($pathways[0]);
		my @literature = ();
		@literature    = (
		    map {   
			s/^\s+//;  # strip leading spaces
			s/\s+$//;  # strip trailing spaces
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
		
		if (defined $entries{$gene}{$to}{ID}){
		    $entries{$gene}{$to}{ID} .= "\n".$goi || '';
		}
		else{
		    $entries{$gene}{$to}{ID} = $goi || 'NA';
		}

		push @{$entries{$gene}{$to}{SYNONYMS}},  @synonyms;		
		@{$entries{$gene}{$to}{SYNONYMS}} = @{unique_array(\@{$entries{$gene}{$to}{SYNONYMS}})};
		push @{$entries{$gene}{$to}{PATHWAY}}, @pathways;
		@{$entries{$gene}{$to}{PATHWAY}} = @{unique_array(\@{$entries{$gene}{$to}{PATHWAY}})};
		push @{$entries{$gene}{$to}{LITERATURE}}, @literature;   
		@{$entries{$gene}{$to}{LITERATURE}} = @{unique_array(\@{$entries{$gene}{$to}{LITERATURE}})};
		my $name = $gene;
		$name = $synonyms[0] if ($name =~ /.goi./ && defined $synonyms[0]);
		if (defined $entries{$gene}{$to}{NAME} && $entries{$gene}{$to}{NAME} ne 'NA'){
		    $entries{$gene}{$to}{NAME} .= "\n".$name || '';
		}
		else{
		    $entries{$gene}{$to}{NAME} = $name || 'NA';
		}

		$entries{$gene}{$to}{TEX}  =	"$goi\/$goi\.tex";

		if ($igvs == 1){
		    push @{$entries{$gene}{$to}{IGV}},"$goi\/snapshots/$goi\_igv1.svg";
		}
		elsif ($igvs == 0){
		    push @{$entries{$gene}{$to}{IGV}},"NONE";
		}
		else{
		    for (1..$igvs){
			push @{$entries{$gene}{$to}{IGV}},"$goi\/snapshots/$goi\_igv$_\.svg";
		    }
		}
		if ($ucscs == 1){
		    push @{$entries{$gene}{$to}{UCSC}},"$goi\/snapshots/$goi\_ucsc1.eps";
		}
		elsif ($ucscs == 0){
		    push @{$entries{$gene}{$to}{UCSC}},"NONE";
		}
		else{
		    for (1..$ucscs){
			push @{$entries{$gene}{$to}{UCSC}},"$goi\/snapshots/$goi\_ucsc$_\.eps";
		    }
		}
		if ($sashimi == 1){
		    push @{$entries{$gene}{$to}{SASHIMI}},"$goi\/snapshots/$goi\_sashimi1.svg";
		}
		elsif ($sashimi == 0){
		    push @{$entries{$gene}{$to}{SASHIMI}},"NONE";
		}
		else{
		    for (1..$sashimi){
			push @{$entries{$gene}{$to}{SASHIMI}},"$goi\/snapshots/$goi\_sashimi$_\.svg";
		    }
		}
		$entries{$gene}{$to}{NOTES} .= $notes || 'NA';
		if (defined $entries{$gene}{$to}{EXTRA} && $entries{$gene}{$to}{EXTRA} ne 'NA'){
		    $entries{$gene}{$to}{EXTRA} .= $extra || '';
		}
		else{
		    $entries{$gene}{$to}{EXTRA} = $extra || 'NA';
		}
		foreach my $syn (@synonyms){
		    next if ($syn eq $gene);
		    $entries{$syn}{$to} = $entries{$gene}{$to};		    
		}   
	    }
	    else{
		foreach my $syn (@synonyms){
		    next if ($syn eq $gene);
		    $entries{$gene}{$to} = $entries{$syn}{$to} if (defined $entries{$syn}{$to} && !defined $entries{$gene}{$to});
		}
		push @process, $line unless (defined $entries{$gene}{$to});
	    }
	}
    }
    elsif ($filetoparse =~ /apg/i){
	print STDERR "Processing APG List!\n";
	while (@process) {
	    my $item = shift(@process);
	    next unless ($item =~ /.apg./);
	    (my $line	  = $item) =~ s/,w+//g;
	    my $to = 'APG';
	    my @fields		  = split(/\;/,$line);
#Folder_ID;Hacker-Name;Gen;Duplicates?;Synonyms;Pathway;Literature;Screen IGV;Scale;Diff. Mock;Diff. Ebola;Diff. Marburg;Profile change;found?;Diff. Mock;Diff. Ebola;Diff. Marburg;Profile Change;Screen UCSC;UCSC Conservation;SNPs;SNPs;Sashimi plot;Intron transcripts;Intron transcripts;Upstream, 5'UTR;Downstream, 3'UTR;Extra Screens and why;Notes;Empty Entries?;Text Draft DONE;done.sh;Questions?;Done Lit Team / Done Questions;FRANZI FINAL CHECK;VERENA;STATE
#hg19.apg.00209;Stefanie W.;ABCA1;0;;;;4;0-376;ed;u17d;dd;0;1;00;00;00;low expression;1;11111;0;0;4;0;0;0;0;0;COOL:  significantly downregulated in EBOV23h (327-390-22);#VALUE!;1;1;0;;1;transport molecules across extra and intracellular membranes;CHECKED
	    my $apg		  = $fields[0];
	    my $hacker		  = $fields[1];
	    next if ($hacker eq 'OPTIONAL' || $hacker eq '');
	    (my $gene		  = $fields[2]) =~ s/^\s+//g;
	    $gene = $apg unless ($gene);
	    next if (defined $entries{$gene}{APG}{ID});
#	    $to = 'GOI' if (exists ($entries{$gene}{GOI}));
	    my $duplicate = $fields[3];   
	    my @synonyms  = ();
	    @synonyms	  = (
		map {   
		    s/^\s+//;  # strip leading spaces
		    s/\s+$//;  # strip trailing spaces
		    $_         # return the modified string
		}
		split(/[,]+|[\s]{2,}|\t/,$fields[4])
		) if ($fields[4] ne '');
	    push @synonyms, $gene unless ($gene =~ /.apg./ || $gene eq '');
	    @synonyms = grep /\S/, @synonyms; ## get rid of empty entries
	    if ($duplicate eq '' || $duplicate == 0){
		my @pathways = ();
		@pathways    = (
		    map {   
			s/^\s+//;  # strip leading spaces
			s/\s+$//;  # strip trailing spaces
			$_         # return the modified string
		    }
		    split(/[,]+|[\s]{2,}|\t/,$fields[5])
		    ) if ($fields[5] ne '');
		push @pathways, 'Unknown' unless ($pathways[0]);
		my @literature = ();
		@literature    = (
		    map {   
			s/^\s+//;  # strip leading spaces
			s/\s+$//;  # strip trailing spaces
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
		
		if (defined $entries{$gene}{$to}{ID}){
		    $entries{$gene}{$to}{ID} .= "\n".$apg || '';
		}
		else{
		    $entries{$gene}{$to}{ID} = $apg || 'NA';
		}

		push @{$entries{$gene}{$to}{SYNONYMS}},  @synonyms;		
		@{$entries{$gene}{$to}{SYNONYMS}} = @{unique_array(\@{$entries{$gene}{$to}{SYNONYMS}})};
		push @{$entries{$gene}{$to}{PATHWAY}}, @pathways;
		@{$entries{$gene}{$to}{PATHWAY}} = @{unique_array(\@{$entries{$gene}{$to}{PATHWAY}})};
		push @{$entries{$gene}{$to}{LITERATURE}}, @literature;   
		@{$entries{$gene}{$to}{LITERATURE}} = @{unique_array(\@{$entries{$gene}{$to}{LITERATURE}})};
		my $name = $gene;
		$name = $synonyms[0] if ($name =~ /.apg./ && defined $synonyms[0]);
		if (defined $entries{$gene}{$to}{NAME} && $entries{$gene}{$to}{NAME} ne 'NA'){
		    $entries{$gene}{$to}{NAME} .= "\n".$name || '';
		}
		else{
		    $entries{$gene}{$to}{NAME} = $name || 'NA';
		}
		$entries{$gene}{$to}{TEX}  =	"$apg\/$apg\.tex";
		if ($igvs == 1){
		    push @{$entries{$gene}{$to}{IGV}},"$apg\/snapshots/$apg\_igv1.svg";
		}
		elsif ($igvs == 0){
		    push @{$entries{$gene}{$to}{IGV}},"NONE";
		}
		else{
		    for (1..$igvs){
			push @{$entries{$gene}{$to}{IGV}},"$apg\/snapshots/$apg\_igv$_\.svg";
		    }
		}
		if ($ucscs == 1){
		    push @{$entries{$gene}{$to}{UCSC}},"$apg\/snapshots/$apg\_ucsc1.eps";
		}
		elsif ($ucscs == 0){
		    push @{$entries{$gene}{$to}{UCSC}},"NONE";
		}
		else{
		    for (1..$ucscs){
			push @{$entries{$gene}{$to}{UCSC}},"$apg\/snapshots/$apg\_ucsc$_\.eps";
		    }
		}
		if ($sashimi == 1){
		    push @{$entries{$gene}{$to}{SASHIMI}},"$apg\/snapshots/$apg\_sashimi1.svg";
		}
		elsif ($sashimi == 0){
		    push @{$entries{$gene}{$to}{SASHIMI}},"NONE";
		}
		else{
		    for (1..$sashimi){
			push @{$entries{$gene}{$to}{SASHIMI}},"$apg\/snapshots/$apg\_sashimi$_\.svg";
		    }
		}
		$entries{$gene}{$to}{NOTES} .= $notes || 'NA';
		if (defined $entries{$gene}{$to}{EXTRA}){
		    $entries{$gene}{$to}{EXTRA} = '' if ($entries{$gene}{$to}{EXTRA} eq 'NA');
		    $entries{$gene}{$to}{EXTRA} .= $extra || '';
		}
		else{
		    $entries{$gene}{$to}{EXTRA} = $extra || 'NA';
		}
		
		foreach my $syn (@synonyms){
		    next if ($syn eq $gene);
		    $entries{$syn}{$to}=$entries{$gene}{$to};		    
		}   
	    }
	    else{
		foreach my $syn (@synonyms){
		    next if ($syn eq $gene);
		    $entries{$gene}{$to} = $entries{$syn}{$to} if (defined $entries{$syn}{$to} && !defined $entries{$gene}{$to});
		}
		push @process, $line unless (defined $entries{$gene}{$to});
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

sub goi_link{
    my $synonym = shift;
    my $goi = shift;
    my $goilink = $goi.".html";
    my $goi_link_html = "<a href=\"$goilink\">$synonym</a>";
    return $goi_link_html;
}
sub index_entry_detailed{
    my $template_path = shift;
    my $name	      = shift;
    my $synonyms      = shift;
    my $goiid	      = shift;
    my $textxt	      = shift;
    my $igv	      = shift;
    my $sashimi	      = shift;
    my $ucsc	      = shift;
    my $maxy	      = shift;
    my ($index_entry,$c3,$c7,$c23,$e3,$e7,$e23,$v3,$v7,$v23);
    my $template     =	Template->new({
	INCLUDE_PATH => ["$template_path"],
	RELATIVE     =>	0
                                      });
    if($maxy ne "NA"){
        $maxy =~ s/\n//g;
        my @max_reads = split (",", $maxy);
        $c3 = $max_reads[0];
        $c7 = $max_reads[1];
        $c23 = $max_reads[2];
        $e3 = $max_reads[3];
        $e7 = $max_reads[4];
        $e23 = $max_reads[5];
        $v3 = $max_reads[6];
        $v7 = $max_reads[7];
        $v23 = $max_reads[8];
    }else{
        $c3 = "NA";
        $c7 = "NA";
        $c23 = "NA";
        $e3 = "NA";
        $e7 = "NA";
        $e23 = "NA";
        $v3 = "NA";
        $v7 = "NA";
        $v23 = "NA";
    }
    my $entry_file = "indexentry.html";
    my $entry_vars = 
    {   
        name	 => $name,
        synonyms => $synonyms,
        goiid	 => $goiid,
        c3 => $c3,
        c7 => $c7,
        c23 => $c23,
        e3 => $e3,
        e7 => $e7,
        e23 => $e23,
        v3 => $v3,
        v7 => $v7,
        v23 => $v23,
        textxt	 => $textxt,
        igv	 => $igv,
        sashimi	 => $sashimi,
        ucsc	 => $ucsc
    };
    $template->process($entry_file,$entry_vars,\$index_entry) || die "Template process failed: ", $template->error(), "\n";
    return $index_entry;
}
sub image_entry{
    my $file = shift;
    my $dir = shift;
    my $odir = shift;
    my $image_entry;
    if($file eq 'NONE'){
        $image_entry = "NONE";
    }else{
        my @file = split ("/", $file);
        my $filename = $file[2];
        if($filename =~/.svg/){
            $filename =~ s/.svg/.png/;
            my $snapshotdir = join("/",$wdir,$dir,$file[0],$file[1]);
            my $imagelink = $wdir . "/" . $dir ."/". $file;
            my $thumblink = "./" ."thumbs/" . "$filename";
            `inkscape --file=$imagelink --export-width=150 --export-height=150 --without-gui --export-png=$thumblink` unless (-e $thumblink || !-e $imagelink);
            $image_entry = "<a href=\"$imagelink\"><img src=\"$thumblink\"></a>";
        }
        if($filename =~/.eps/){
            $filename =~ s/.eps/.png/;
            my $snapshotdir = join("/",$wdir,$dir,$file[0],$file[1]);
            my $imagelink = $wdir . "/" . $dir ."/". $file;
            my $thumblink = "./" ."thumbs/" . "$filename";
            `convert $imagelink -resize 150x150! $thumblink` unless (-e $thumblink || !-e $imagelink);
            $image_entry = "<a href=\"$imagelink\"><img src=\"$thumblink\"></a>";
        }
    }
    return $image_entry;
}

sub additional_plot_entry{
    my $additional_plot_entry = shift;
    my $file = shift;
    my $dir = shift;
    my $link_entry;
    if($additional_plot_entry eq '' || $additional_plot_entry eq 'NA' || $additional_plot_entry =~ /^0/ || $file eq 'NONE'){
        $link_entry = "NA";
    }else{
        my @file = split ("/", $file);
	my $snapshotdir = join("/",$wdir,$dir,$file[0],$file[1]);
        $link_entry = "<a href=\"$snapshotdir\">Additional plots</a>";
    }
    return $link_entry;
}

sub link_entry{
    my $file = shift;
    my $dir = shift;
    my $link_entry;
    if($file eq 'NONE'){
        $link_entry = "NONE";
    }else{
        $file =~ s/.tex/.pdf/;
        my $filelink = $wdir . "/" . $dir ."/". $file;
        $link_entry = "<a href=\"$filelink\">pdf</a>";
    }
    return $link_entry;
}

