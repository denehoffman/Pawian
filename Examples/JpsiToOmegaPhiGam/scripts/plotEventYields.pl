#!/usr/bin/perl

my $baseDir = "/data/sleipnir1/marc/pascal/MassIndepFits/";
my $criterion = @ARGV[0]; #"AICc";
print $criterion ."\n";

my $minMass=1.8;
my $binWidth=0.05;
my $nbins=25;

my %allWaves;
my @massBinResults;

for(my $bin=0;$bin<$nbins;$bin++){
    
    $mass=($minMass+($bin+0.5)*$binWidth)*1000;
    my $bestHypo="";
    
    my $listByCriterion=$criterion."List" .$mass .".txt";
    
    open (LIST, $listByCriterion);
    while(<LIST>){
	($bestHypo) = split;
	last;
    }
    close LIST;
    
    die("Failed to extract best hypothesis from $listByCriterion\n") if($bestHypo eq "");
    
    print "The best hypotheis is $bestHypo\n";
    
    my %waves;
    extractWaveYields($bestHypo, $mass, \%waves);
    print "Final\n";
    foreach(keys(%waves)){
	$allWaves{$_}=1;
    }
    push @massBinResults, \%waves;
}

my %colors =("F1"=>"kBlue", 
	    "Eta"=>"kRed" ,
	    "F0"=>"kGreen",
	    "F2"=>"kBlack",
	    "Eta2"=>"kCyan" );


open PLOT, ">plot_".$criterion.".C";
print PLOT "void plot_".$criterion."(){\n";

foreach(keys %allWaves){
    my $histo="h_$_";
    my $maxMass=$minMass+$nbins*$binWidth;
    print PLOT "TH1F* $histo = new TH1F(\"$histo\",\"$histo\",$nbins,$minMass,$maxMass);\n";
    print PLOT "$histo->SetMinimum(0);\n";
    print PLOT "$histo->SetMarkerStyle(20);\n";
    print PLOT "$histo->SetMarkerColor($colors{$_});\n";
}

$bin=1;
foreach my $ref(@massBinResults){
    foreach my$wave(keys %allWaves){
	my $histo="h_$wave";

	if($ref->{$wave}){
	    print PLOT $histo ."->SetBinContent($bin, $ref->{$wave});\n"; 
	    print PLOT $histo ."->SetBinError($bin, 1);\n"; 
	}else{
	    print PLOT $histo ."->SetBinContent($bin, 0.0);\n";
	    print PLOT $histo ."->SetBinError($bin, 1);\n"; 
	}
    }
    $bin++;
}

my $first=1;
foreach(keys %allWaves){
    my $histo="h_$_";
    if($first){
	print PLOT $histo ."->GetXaxis()->SetRangeUser(0,900);\n";
	print PLOT $histo ."->Draw(\"E\");\n";
	$first=0;
	next;
    }
    print PLOT $histo ."->Draw(\"Esame\");\n";
    
}

print PLOT "}\n";

close PLOT;


sub extractWaveYields(){
    my $bestHypo=shift;
    my $mass =shift;
    my $waveRef = shift;
    
    foreach my $theWave(split "_", $bestHypo){
	my $waveYield=0;
	my $theQaLogFile = "$baseDir/$mass/$bestHypo/qa/$theWave/qa.log";
	print $theQaLogFile."\n";
	open(LOGFILE,$theQaLogFile);
	while(<LOGFILE>){
	    if(/No of fit events\s+(\d.+\d)/){
		$waveYield=$1;
		last;
	    }
	}
	close LOGFILE;
	$waveRef->{$theWave}=$waveYield;
    }
}




