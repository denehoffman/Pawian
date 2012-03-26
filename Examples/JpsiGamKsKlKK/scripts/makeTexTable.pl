#!/usr/bin/perl -w

use strict;

sub sortParams {
    my $a = $_[0];
    my $b = $_[1];
    
    my %dlsa = decayLS($a);
    my %dlsb = decayLS($b);
    
    if($dlsa{"state"} eq $dlsb{"state"} && $dlsa{"prod"} eq $dlsb{"prod"}){
	if($dlsa{"S"} <=> $dlsb{"S"} ){ return ($dlsa{"S"} <=> $dlsb{"S"}) };
	if($dlsa{"L"} <=> $dlsb{"L"} ){ return ($dlsa{"L"} <=> $dlsb{"L"}) };
    }
    if( $dlsa{"state"} eq $dlsb{"state"} && $dlsa{"prod"} <=> $dlsb{"prod"}  ) {return -($dlsa{"prod"} <=> $dlsb{"prod"})};
    return ($dlsa{"state"} cmp $dlsb{"state"} );
}


sub decayLS {
    my $pname = $_[0];
    
    my %dls;
    if($pname=~/L(.)S(.)/){
	$dls{"L"}=$1;
	$dls{"S"}=$2;
    }
    if($pname=~/PsiTo(.+)Gam/){
	$dls{"state"} = $1;
	$dls{"prod"} = 1;
    }
    if($pname=~/S.(.+)ToPhiPhi/){
	$dls{"state"} = $1;
	$dls{"prod"} = 0;
    }
    return %dls;
}



my %params;
my $resultFile = $ARGV[0];
open(FILE, $resultFile);
foreach my $line(<FILE>){
    chomp $line;
    if($line=~ /(^J.+)\t(.+)\t(.+)/){
	my $parname =$1;
	my $parval = $2;
	my $parerr = $3;
	
	my $pair = [];
	push @$pair, $parval;
	push @$pair, $parerr;
	
	$params{$parname} = $pair;
    }
}
close FILE;

open TEX, ">tmp.tex";
print TEX "\\documentclass{article}\n";
print TEX "\\begin{document}\n";
print TEX "\\begin{verbatim}$resultFile\\end{verbatim}\n\n";
print TEX "\\begin{table}\n";
print TEX "\\begin{tabular}{l|l|l|l|l}\n";
print TEX "Ampl. & M & dM & Phi & dPhi \\\\ \n";
foreach my $parname( sort {sortParams($a,$b)} (keys %params) ){
    
    next if( !($parname=~/Mag/) );
    my $magName = $parname;
    my $phiName = $parname;
    my $name = $parname;
    
    my %dls = decayLS($parname);
    
    $name=~s/Mag//;
    $phiName=~s/Mag/Phi/;
    
    my $mag = $params{$magName}[0];
    my $magErr = $params{$magName}[1];
    
    my $phi = $params{$phiName}[0];
    my $phiErr = $params{$phiName}[1];
    
    my $tname = $dls{"state"} ." (L=" .$dls{"L"} ." S=" .$dls{"S"} .") " .$dls{"prod"};
    my $ls = "(L=" .$dls{"L"} ." S=" .$dls{"S"} .")";
    my $pflag="D";
    if($dls{"prod"}){$pflag="P"};
    $tname = sprintf("%6s %s %s",$dls{"state"},$pflag, $ls);
    printf("%s\t%.3f\t%.3f\t%.3f\t%.3f\n", $tname, $mag, $magErr, $phi, $phiErr);
    
    printf TEX ("%s&%.3f&%.3f&%.3f&%.3f\\\\\n", $tname, $mag, $magErr, $phi, $phiErr);
    
}
print TEX "\\end{tabular}\n\\end{table}\n\\end{document}\n";

close TEX;

system("pdflatex tmp.tex; cp tmp.pdf ~/WWW/d0.pdf");

