#!/usr/bin/perl -w

use strict;
use File::Basename;


my $baseDir = $ARGV[0];

my $templConfig = $baseDir."/MassIndep.cfg";
my $templStart = $baseDir."/startParams.dat";

my $queue="pwa";
my $app = "/data/sleipnir1/marc/PawianPhase/Pawian/Examples/JpsiGamKsKlKK/bin/gcc-4.1.2/release/link-static/MJpsiGamKlKsKKApp";




my $newConfig = `cat $templConfig`;

my $massMin=2.05;
my $binWidth=0.1;
my $nbins=10;
my $bin=0;
while($bin<$nbins){
    
    my $rangeMin=$massMin+$bin*$binWidth;
    my $rangeMax = $rangeMin+$binWidth;
    my $rangeMean = sprintf("%.0f",($rangeMin+0.5*$binWidth)*1000);
    
    my $fitDir = $baseDir."/" .$rangeMean;
    my $fitConfig = $fitDir ."/" .basename($templConfig);
    my $fitStart = $fitDir ."/" .basename($templStart);
    
    if(! (-e $fitDir)) {system("mkdir $fitDir");}
    open CONF, ">$fitConfig";
    foreach my $line(split "\n", $newConfig){
	if($line=~/massRangeMin/){
	    print CONF "massRangeMin = $rangeMin\n";
	    next; #line
	}
	if($line=~/massRangeMax/){
            print CONF "massRangeMax = $rangeMax\n";
            next; #line
        }
	if($line=~/paramFile/){
	    print CONF "paramFile = $fitStart\n";
	    next; #line
	}
	print CONF $line ."\n";
    }
    close CONF;
    
    
    system("cp $templStart $fitDir/");
    chdir($fitDir);
    my $cmd = "bsub -q $queue -o $fitDir/batch.log \"$app -c $fitConfig --mode pwa >&! $fitDir/fit.log\"";
    print $cmd."\n";
    system($cmd);
    chdir("-");

    $bin++;
}
