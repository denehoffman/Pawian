#!/usr/bin/perl -w

use strict;
use File::Basename;


my $baseDir = "/data/sleipnir1/marc/pascal/MassIndepFits/";
print $baseDir."\n";

my $pawiandir = $ENV{CMAKE_SOURCE_DIR};
print $pawiandir."\n";


my $queue="batch";

my $app = "$pawiandir/Examples/JpsiToOmegaPhiGam/bin/gcc-4.4.6/release/link-static/MJpsiToOmegaPhiGamApp";
my @mySingleHypList=('All', 'Eta2', 'Eta', 'F0', 'F1', 'F2');
my @myHypList=('Eta', 'Eta2', 'F0', 'F1', 'F2', 'Eta_F0', 'Eta_F1', 'Eta_F2', 'Eta_Eta2', 'F0_F1', 'F0_F2', 'F0_Eta2', 'F1_F2', 'F1_Eta2', 'F2_Eta2', 'Eta_F0_F1', 'Eta_F0_F2', 'Eta_F0_Eta2', 'Eta_F1_F2', 'Eta_F1_Eta2', 'Eta_F2_Eta2', 'F0_F1_F2', 'F0_F1_Eta2', 'F0_F2_Eta2', 'F1_F2_Eta2', 'Eta_F0_F1_F2', 'Eta_F0_F1_Eta2', 'Eta_F0_F2_Eta2', 'Eta_F1_F2_Eta2', 'F0_F1_F2_Eta2', 'Eta_F0_F1_F2_Eta2'  ); 


my $theHyp= "";

for $theHyp (@myHypList) {
    my $configFile = "JpsiGamOmegaPhi_$theHyp.cfg";

my $allHyps = "All";

my $massMin=1.8;
my $binWidth=0.05;
my $nbins=25;
my $bin=0;

while($bin<$nbins){
    
    my $rangeMin=$massMin+$bin*$binWidth;
    my $rangeMax = $rangeMin+$binWidth;
    my $rangeMean = sprintf("%.0f",($rangeMin+0.5*$binWidth)*1000);
     
    my $rangeDir= $baseDir."/" .$rangeMean;

    if(! (-e $rangeDir)){
	print "rangeDir:\t".$rangeDir."\tdoes not exist!\n";
	exit(0);
    } 
 
    my $fitDir = $rangeDir."/" .$theHyp;
    
    if(! (-e $fitDir)) {
	print "fitDir:\t".$fitDir."\tdoes not exist!\n";
	exit(0);
    }

my $fitConfig = $fitDir."/JpsiGamOmegaPhi_$theHyp.cfg";
print "fitConfig:\t".$fitConfig."\n"; 

my $qaDir = $fitDir."/qa";
if(! (-e $qaDir)) {system("mkdir $qaDir");}

my $singleHyp= "";
    my $theCompHyps = $theHyp;

for $singleHyp (@mySingleHypList) {
    if ($theCompHyps =~ /$singleHyp/ || $singleHyp =~ /All/){
	if ($singleHyp !~ /All/){ $theCompHyps =~ s/$singleHyp/''/g};
	    print $singleHyp."\t exists\n";
	my $qaHypDir = $qaDir."/".$singleHyp;
	if(! (-e $qaHypDir)) {system("mkdir $qaHypDir");}
#	chdir($qaHypDir);
        my $qaConfigFile =  $qaHypDir."/JpsiGamOmegaPhi_$theHyp.cfg";
#	system("cp $fitConfig $qaConfigFile");

	open INF, "<$fitConfig";
	open OUTF, ">$qaConfigFile";
	while (<INF>) {
            if ($_ !~ /paramFile/ && $_ !~ /enableHyp/) { print OUTF $_; }
            elsif ($_ =~ /paramFile/) { print OUTF "paramFile = ".$fitDir."/finalResult.dat"; }
	    elsif($_ =~ /enableHyp/){
		my $tmpRegexp = $singleHyp."_";
		if ($singleHyp =~ /All/ || $_ =~ /$tmpRegexp/) { print OUTF $_; }
	    }
	}
	close INF;
	close OUTF;
	
	chdir($qaHypDir);
	next if(-e "$qaHypDir/qaSummary.dat");
	my $cmd = "bsub -q $queue -o $qaHypDir/batch.log \"$app -c $qaConfigFile  --mode qaMode >&! $qaHypDir/qa.log\"";
	print $cmd."\n";
	system($cmd);
    }
}

    chdir($baseDir);

    $bin++;
}
}
