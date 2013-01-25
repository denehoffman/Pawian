#!/usr/bin/perl -w

use strict;
use File::Basename;

open FIXPARAMFILE, "FixParams.txt" or die $!;

my @FixParamMap = map { chomp; [split /\t/]; } <FIXPARAMFILE>;

# print the @AOA array
foreach my $item1 (@FixParamMap){
    print "@{$item1}[0]\n";
    print "@{$item1}[1]\n";
    print "@{$item1}[2]\n";
  print "\n";
}

#exit(0);

my $baseDir = "/data/sleipnir1/marc/pascal/MassIndepFits/";
print $baseDir."\n";

my $pawiandir = $ENV{CMAKE_SOURCE_DIR};
print $pawiandir."\n";

my $configFilesDir = $pawiandir."/Examples/JpsiToOmegaPhiGam/ConfigFiles";
print $configFilesDir."\n"; 

#my $queue="users";
my $app = "$pawiandir/Examples/JpsiToOmegaPhiGam/bin/gcc-4.4.6/release/link-static/MJpsiToOmegaPhiGamApp";

my @myHypList=('Eta', 'Eta2', 'F0', 'F1', 'F2', 'Eta_F0', 'Eta_F1', 'Eta_F2', 'Eta_Eta2', 'F0_F1', 'F0_F2', 'F0_Eta2', 'F1_F2', 'F1_Eta2', 'F2_Eta2', 'Eta_F0_F1', 'Eta_F0_F2', 'Eta_F0_Eta2', 'Eta_F1_F2', 'Eta_F1_Eta2', 'Eta_F2_Eta2', 'F0_F1_F2', 'F0_F1_Eta2', 'F0_F2_Eta2', 'F1_F2_Eta2', 'Eta_F0_F1_F2', 'Eta_F0_F1_Eta2', 'Eta_F0_F2_Eta2', 'Eta_F1_F2_Eta2', 'F0_F1_F2_Eta2', 'Eta_F0_F1_F2_Eta2'  );


my $theHyp= "";

for $theHyp (@myHypList) {
    print "Hyp:\t$theHyp\n";
    my $addParFixList = "";
    my $replaceHypList ='';
    my @currentHyps  = split(/_/,$theHyp);
    foreach my $currentHyp (@currentHyps){
#	print "currentHyp:\t$currentHyp\n";
        my $currentHypComp=$currentHyp."_";
	foreach my $item1 (@FixParamMap){
	    my $CurrenParamName = @{$item1}[0];
#	    print "CurrenParamName:\t$CurrenParamName\n";
            if (index($CurrenParamName, $currentHypComp) != -1){
#		print "$CurrenParamName\n";
#		print "@{$item1}[1]\n";
#		print "@{$item1}[2]\n";
#		print "\n";
		$addParFixList ="$addParFixList\nmnParFix\t=\t$CurrenParamName";
		$replaceHypList = "$replaceHypList\n$CurrenParamName";
	    }
	}
    }

    my $firstHyp= $currentHyps[0];
    print "first Hyp:\t".$firstHyp."\n";
print "***Fix List***\n";
print "$addParFixList\n";
    print "***replace list***\n";
    print "$replaceHypList\n";
#next;

my $templConfig = $configFilesDir."/JpsiToOmegaPhiGamDefault.cfg";
my $templStart = $configFilesDir."/startParams$theHyp.cfg";

my $newConfig = `cat $templConfig`;
my $newConfigStart = `cat $templStart`;

my $massMin=1.8;
my $binWidth=0.05;
my $nbins=25;
my $bin=0;
while($bin<$nbins){
    
    my $rangeMin=$massMin+$bin*$binWidth;
    my $rangeMax = $rangeMin+$binWidth;
    my $rangeMean = sprintf("%.0f",($rangeMin+0.5*$binWidth)*1000);
     
    my $rangeDir= $baseDir."/" .$rangeMean;
    if(! (-e $rangeDir)) {system("mkdir $rangeDir");}
 
    my $fitDir = $rangeDir."/" .$theHyp;
    my $fitConfig = $fitDir ."/JpsiGamOmegaPhi_$theHyp.cfg";
    my $fitStart = $fitDir."/" .basename($templStart);  
    
    my $finalResult =  $fitDir ."/finalResult.dat";
    if( -e $finalResult){
	print "Fit $finalResult allready finished successfully -- do not submit again!\n";
	$bin++;
	next;
    }
    

    
    if(! (-e $fitDir)) {system("mkdir $fitDir");}
    open CONF, ">$fitConfig";
#    print CONF $addParFixList."\n";
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
    print CONF "mnParFix = J1P-1C-1Lama0Lamb1Gamma".$firstHyp."_2000Phi\n";
    print CONF $addParFixList."\n";
    foreach my $currentHyp (@currentHyps){
#	print "currentHyp:\t$currentHyp\n";
        my $currentHypComp=$currentHyp."_";
	print CONF "enableHyp = Gamma".$currentHyp."_2000\n";
    }
    close CONF;
    
    open CONFSTART, ">$fitStart";
    foreach my $line(split "\n", $newConfigStart){
        my @lineCont = split(/\t/, $line);
        my $hypAtLine=$lineCont[0];
        my $replaceLine=0;
        my $replaceLineCont = "";
#	foreach my $currentReplaceHyp(split "\n", $replaceHypList){
	foreach my $currentReplaceHypLine (@FixParamMap){
#	    my @currentReplaceHypLineCont = @{$currentReplaceHypLine}[0];   
	    my $currentReplaceHyp = @{$currentReplaceHypLine}[0];
#            print "currentReplaceHyp: ".$currentReplaceHyp."\n";
	    if ($currentReplaceHyp =~ $hypAtLine){
		$replaceLine=1;
		$replaceLineCont = @{$currentReplaceHypLine}[0]."\t". @{$currentReplaceHypLine}[1]."\t". @{$currentReplaceHypLine}[2]; 
	    } 
	}
	if ($replaceLine==0) {print CONFSTART $line ."\n";}
	else {print CONFSTART $replaceLineCont."\n"}
    }
    close CONFSTART;

#    system("cp $templStart $fitDir/");
    chdir($fitDir);
    my $cmd = "bsub -l nodes=1:ppn=16 -q pwa -o $fitDir/batch.log \"$app -c $fitConfig  --mode pwa >&! $fitDir/fit.log\"";
    #my $cmd = "bsub -l nodes=1:ppn=4 -q users  -o $fitDir/batch.log \"$app -c $fitConfig  --mode pwa >&! $fitDir/fit.log\"";
    print $cmd."\n";
    system($cmd);
#    chdir("-");
    chdir($baseDir);

    $bin++;
}
}
