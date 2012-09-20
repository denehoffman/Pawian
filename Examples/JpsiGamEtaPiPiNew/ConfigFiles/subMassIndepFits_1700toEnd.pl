#!/usr/bin/perl -w

use strict;
use File::Basename;

my $baseDir = $ARGV[0];
print $baseDir."\n";

my $pawiandir = $ENV{CMAKE_SOURCE_DIR};
print $pawiandir."\n";

my $configFilesDir = $pawiandir."/Examples/JpsiGamEtaPiPiNew/ConfigFiles";
my $summaryFileDir = $pawiandir."/Examples/JpsiGamEtaPiPiNew/bin/gcc-4.4.4/release/link-static/scan/summary.dat";
my $mcFileDir = "/data/liema/michaell/PwaJpsiGamEtaPiPi/FS4Vectors-PHSPMC-JpsiGamEtaPipPim_4m_v4.dat";

print $configFilesDir."\n"; 

my $queue="pwa";
my $app = "$pawiandir/Examples/JpsiGamEtaPiPiNew/bin/gcc-4.4.4/release/link-static/MJpsiGamEtaPiPiNewApp";

#my @myHypList=('Eta', 'Eta2', 'F1', 'EtaF1', 'EtaEta2', 'EtaF1Eta2');
my @myHypList=('Eta2');

#my @mySubDecayList=('_a0f0f2sigma', '_a0f0a2f2', '_a0f0a2sigma', '_a0f0a2f2sigma');
my @mySubDecayList=('_a0f0a2f2', '_a0f0a2sigma', '_a0f0a2f2sigma');

my $theHyp= "";
my $subDecay="";

for $theHyp (@myHypList) {
    for $subDecay (@mySubDecayList){
	my $templConfig = $configFilesDir."/JpsiGamEtaPiPiNew_$theHyp$subDecay.cfg";
	my $templStart = $configFilesDir."/defaultparams.dat";
	
	my $newConfig = `cat $templConfig`;
	
	my $massMin=1.70;
	my $binWidth=0.04;
	my $nbins=8;
	my $bin=0;
	while($bin<$nbins){
	    
	    my $rangeMin=$massMin+$bin*$binWidth;
	    my $rangeMax = $rangeMin+$binWidth;
	    my $rangeMean = sprintf("%.0f",($rangeMin+0.5*$binWidth)*1000);
	
	    my $rangeDir= $baseDir."/" .$rangeMean;
	    if(! (-e $rangeDir)) {system("mkdir $rangeDir");}
	    
	    my $fitDir = $rangeDir."/" .$theHyp;
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
		if($line=~/sumFile/){
                    print CONF "sumFile = $summaryFileDir\n";
                    next; #line
                }
		if($line=~/mcFile/){
                    print CONF "mcFile = $mcFileDir\n";
                    next; #line
                }
		
		print CONF $line ."\n";
	    }
	    close CONF;
	    
	    
	    system("cp $templStart $fitDir/");
	    chdir($fitDir);
	    my $cmd = "bsub -l nodes=1:ppn=16 -q $queue -o $fitDir/batch.log \"$app -c $fitConfig  --mode pwa >&! $fitDir/log$subDecay.out \"";
	    print $cmd."\n";
	    system($cmd);
#    chdir("-");
	    chdir($baseDir);
	    
	    $bin++;
	}
    }

}
