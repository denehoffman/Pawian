#!/usr/bin/perl -w

use strict;
use File::Basename;



my $baseDir = "/data/sleipnir1/marc/pascal/MassIndepFits";

my @myHypList=('Eta', 'Eta2', 'F0', 'F1', 'F2', 'Eta_F0', 'Eta_F1', 'Eta_F2', 'Eta_Eta2', 'F0_F1', 'F0_F2', 'F0_Eta2', 'F1_F2', 'F1_Eta2', 'F2_Eta2', 'Eta_F0_F1', 'Eta_F0_F2', 'Eta_F0_Eta2', 'Eta_F1_F2', 'Eta_F1_Eta2', 'Eta_F2_Eta2', 'F0_F1_F2', 'F0_F1_Eta2', 'F0_F2_Eta2', 'F1_F2_Eta2', 'Eta_F0_F1_F2', 'Eta_F0_F1_Eta2', 'Eta_F0_F2_Eta2', 'Eta_F1_F2_Eta2', 'F0_F1_F2_Eta2', 'Eta_F0_F1_F2_Eta2');
 

my $qaFile = "qaSummary.dat";

my $massMin=1.80;
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
    
    my $summaryFile = $baseDir."/summary".$rangeMean.".txt";
    open SUMFILE, ">$summaryFile";
    print SUMFILE "Phi Phi Mass range\tMin: ".$rangeMin."\tMax: ".$rangeMax."\n\n";
#    print SUMFILE "hypothesis\tBIC\tAIC\tNLL\tNoParams\n\n";
    printf SUMFILE '%17s',"hypothesis";
    printf SUMFILE '%10s',"BIC";
    printf SUMFILE '%10s',"AIC";
    printf SUMFILE '%10s',"AICc";
    printf SUMFILE '%10s',"NLL";
    printf SUMFILE '%12s',"NoParams\n";
 
    my $makeSummaryListFile = $baseDir."/makeSummaryLists".$rangeMean.".C";
    open MAKESUMMARYLISTFILE, ">$makeSummaryListFile";

    print MAKESUMMARYLISTFILE "#include <iostream>\n";
    print MAKESUMMARYLISTFILE "#include <fstream>\n";
    print MAKESUMMARYLISTFILE "#include <map>\n";
    print MAKESUMMARYLISTFILE "#include <string>\n";
    print MAKESUMMARYLISTFILE "#include \"TMath.h\"\n";
    print MAKESUMMARYLISTFILE "#include \"Math/DistFunc.h\"\n";

    print MAKESUMMARYLISTFILE "using namespace std;\n";

    print MAKESUMMARYLISTFILE "void makeSummaryLists$rangeMean\(\)\{\n\n";
    print MAKESUMMARYLISTFILE "std::map< string, pair<double, int> > lhResultMap;\n";
    print MAKESUMMARYLISTFILE "std::map< double, string > bicMap;\n";
    print MAKESUMMARYLISTFILE "std::map< double, string > aicMap;\n\n";
    print MAKESUMMARYLISTFILE "std::map< double, string > aiccMap;\n\n";
    my $theHyp= "";
    for $theHyp (@myHypList) {    

    my $fitDir = $rangeDir."/" .$theHyp;
    
    if(! (-e $fitDir)) {
        print "fitDir:\t".$fitDir."\tdoes not exist!\n";
        exit(0);
    }

    my $qaDir = $fitDir."/qa/All";
    if(! (-e $qaDir)) {
        print "qaDir:\t".$qaDir."\tdoes not exist!\n";
        exit(0);
    }

    my $qaFile = $qaDir."/qaSummary.dat";
    my $qaCont = `cat $qaFile`;
    
    foreach my $line(split "\n", $qaCont){
	if($line=~/BIC/){
	    my @lineCont = split(/\t/,$line);
#           print SUMFILE  $theHyp."\t".$lineCont[1]."\t";
	    printf SUMFILE '%17s',"$theHyp";
	    printf SUMFILE '%10s',"$lineCont[1]";
            printf MAKESUMMARYLISTFILE "bicMap\[$lineCont[1]\]=\"$theHyp\";\n";
	}
	elsif($line=~/AICa/){
	    my @lineCont = split(/\t/,$line);
#            print SUMFILE  $lineCont[1]."\t";
            printf SUMFILE  '%10s',"$lineCont[1]";
	    printf MAKESUMMARYLISTFILE "aicMap\[$lineCont[1]\]=\"$theHyp\";\n";
	}
	elsif($line=~/AICc/){
	    my @lineCont = split(/\t/,$line);
#            print SUMFILE  $lineCont[1]."\t";
            printf SUMFILE  '%10s',"$lineCont[1]";
	    printf MAKESUMMARYLISTFILE "aiccMap\[$lineCont[1]\]=\"$theHyp\";\n";
	}
	elsif($line=~/logLh/){
	    my @lineCont = split(/\t/,$line);
#           print SUMFILE  $lineCont[1]."\t";
            printf SUMFILE  '%10s',"$lineCont[1]"; 
	    printf MAKESUMMARYLISTFILE "lhResultMap[\"$theHyp\"\] = make_pair($lineCont[1],";
	}
	elsif($line=~/free parameter/){
	    my @lineCont = split(/\t/,$line);
#            print SUMFILE  $lineCont[1]."\n";
            printf SUMFILE  '%12s',"$lineCont[1]\n";
	    printf MAKESUMMARYLISTFILE "$lineCont[1]\);\n\n";
	}
    }
  }

  print MAKESUMMARYLISTFILE "std::map<double, string> ratioResultMapInv;\n";  
  print MAKESUMMARYLISTFILE "std::map< string, pair<double, int> >::iterator it;\n";
  print MAKESUMMARYLISTFILE "for (it=lhResultMap.begin(); it!=lhResultMap.end(); ++it){\n";
  print MAKESUMMARYLISTFILE "double currentLhRatio=sqrt(ROOT::Math::chisquared_quantile_c(TMath::Prob(2.*fabs(lhResultMap[\"Eta_F0_F1_F2_Eta2\"].first-it->second.first), lhResultMap[\"Eta_F0_F1_F2_Eta2\"].second-it->second.second),1));\n";
    print MAKESUMMARYLISTFILE     "cout << \"LhRatio:\t\" << it->first << \"\t\" << currentLhRatio << endl;\n";
    print MAKESUMMARYLISTFILE  "ratioResultMapInv[currentLhRatio]=it->first;\n";
    print MAKESUMMARYLISTFILE "\}\n\n";
    print MAKESUMMARYLISTFILE "ofstream lhStream;\n";
    print MAKESUMMARYLISTFILE "lhStream.open(\"LhRatioList$rangeMean.txt\");\n";
    print MAKESUMMARYLISTFILE "map<double, string>::iterator it2;\n";
    print MAKESUMMARYLISTFILE "for (it2=ratioResultMapInv.begin(); it2!=ratioResultMapInv.end(); ++it2){\n";
    print MAKESUMMARYLISTFILE "cout << it2->second << \"\t\" << it2->first << endl;\n";
    print MAKESUMMARYLISTFILE "lhStream << it2->second << \"\t\" << it2->first << endl;\n";
    print MAKESUMMARYLISTFILE "}\n\n";

    print MAKESUMMARYLISTFILE "ofstream aicStream;\n";
    print MAKESUMMARYLISTFILE "aicStream.open(\"AICList$rangeMean.txt\");\n";
    print MAKESUMMARYLISTFILE "map<double, string>::iterator it3;\n";
    print MAKESUMMARYLISTFILE "for (it3=aicMap.begin(); it3!=aicMap.end(); ++it3){\n";
    print MAKESUMMARYLISTFILE "cout << it3->second << \"\t\" << it3->first << endl;\n";
    print MAKESUMMARYLISTFILE "aicStream << it3->second << \"\t\" << it3->first << endl;\n";
    print MAKESUMMARYLISTFILE "}\n\n";

    print MAKESUMMARYLISTFILE "ofstream aiccStream;\n";
    print MAKESUMMARYLISTFILE "aiccStream.open(\"AICcList$rangeMean.txt\");\n";
    print MAKESUMMARYLISTFILE "for (it3=aiccMap.begin(); it3!=aiccMap.end(); ++it3){\n";
    print MAKESUMMARYLISTFILE "cout << it3->second << \"\t\" << it3->first << endl;\n";
    print MAKESUMMARYLISTFILE "aiccStream << it3->second << \"\t\" << it3->first << endl;\n";
    print MAKESUMMARYLISTFILE "}\n\n";

    print MAKESUMMARYLISTFILE "ofstream bicStream;\n";
    print MAKESUMMARYLISTFILE "bicStream.open(\"BICList$rangeMean.txt\");\n";
    print MAKESUMMARYLISTFILE "for (it3=bicMap.begin(); it3!=bicMap.end(); ++it3){\n";
    print MAKESUMMARYLISTFILE "cout << it3->second << \"\t\" << it3->first << endl;\n";
    print MAKESUMMARYLISTFILE "bicStream << it3->second << \"\t\" << it3->first << endl;\n";
    print MAKESUMMARYLISTFILE "}\n\n";

    print MAKESUMMARYLISTFILE "}\n";
    close SUMFILE;
    close MAKESUMMARYLISTFILE;
    $bin++;
}

#system("cd $baseDir; find -name \"make*C\" | sed -e \'s/C/C++/\' | xargs -n1 root -b -q -x");

