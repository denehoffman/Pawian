#!/usr/bin/perl -w
use strict;

my $baseDir = "/data/sleipnir1/marc/pascal/MassIndepFits/";
print $baseDir."\n";

my $pawiandir = $ENV{CMAKE_SOURCE_DIR};
print $pawiandir."\n";

my $configFilesDir = $pawiandir."/Examples/JpsiToOmegaPhiGam/ConfigFiles";
print $configFilesDir."\n"; 

my $app = "$pawiandir/Examples/JpsiToOmegaPhiGam/bin/gcc-4.4.6/release/link-static/MJpsiToOmegaPhiGamApp";


my @myHypList=('Eta', 'Eta2', 'F0', 'F1', 'F2', 'Eta_F0', 'Eta_F1', 'Eta_F2', 'Eta_Eta2', 'F0_F1', 'F0_F2', 'F0_Eta2', 'F1_F2', 'F1_Eta2', 'F2_Eta2', 'Eta_F0_F1', 'Eta_F0_F2', 'Eta_F0_Eta2', 'Eta_F1_F2', 'Eta_F1_Eta2', 'Eta_F2_Eta2', 'F0_F1_F2', 'F0_F1_Eta2', 'F0_F2_Eta2', 'F1_F2_Eta2', 'Eta_F0_F1_F2', 'Eta_F0_F1_Eta2', 'Eta_F0_F2_Eta2', 'Eta_F1_F2_Eta2', 'F0_F1_F2_Eta2', 'Eta_F0_F1_F2_Eta2'  );


my $theHyp= "";



for $theHyp (@myHypList) {
    print "Hyp:\t$theHyp\n";
    
    my $tmpConfigFile = "tmpConfigFile.cfg";
    
    my $theConfig = `cat  $configFilesDir/JpsiToOmegaPhiGamDefault.cfg`;
    foreach my $theReso(split "_", $theHyp){
	$theConfig.="enableHyp = Gamma$theReso"."_2000\n";
    }
    open (TMP, ">$tmpConfigFile");
    print TMP 	$theConfig;
    close TMP;
    
    my $paramFile = "$configFilesDir/startParams$theHyp".".cfg";
    
    my $cmd ="$app -c $tmpConfigFile --mode dumpDefaultParams";
    system("$cmd");
    system("mv defaultparams.dat $paramFile");
    
    system("rm $tmpConfigFile");
    
}
