#!/usr/bin/perl

#************************************************************************#
#                                                                        #
#  Copyright 2014 Bertram Kopf (bertram@ep1.rub.de)                      #
#                 Julian Pychy (julian@ep1.rub.de)                       #
#                 - Ruhr-Universität Bochum                              #
#                                                                        #
#  This file is part of Pawian.                                          #
#                                                                        #
#  Pawian is free software: you can redistribute it and/or modify        #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  (at your option) any later version.                                   #
#                                                                        #
#  Pawian is distributed in the hope that it will be useful,             #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with Pawian.  If not, see <http:#www.gnu.org/licenses/>.        #
#                                                                        #
#************************************************************************#

#######################################################
##                                                   ##
##  pawianfixparameters                              ##
##                                                   ##
##  This script generated the mnParFix list to be    ##
##  used in the Pawian configuration file.           ##
##                                                   ##
##  - Julian Pychy     julian@ep1.rub.de             ##
##                                                   ##
#######################################################

use strict;
use warnings;
use Getopt::Long;

my $showHelp=0;
my $filename="";
my $prio="";

GetOptions ("help|h" => \$showHelp,
	    "file=s" => \$filename,
	    "prio=s" => \$prio)  or die "Bad arguments!";

if($showHelp){
    Help();
    exit;
}

if($filename eq "" || !(-e $filename)){
    die "File not found.";
}


# Read the parameter file
open FILE, "<", $filename or die $!;

my @linesInFile;
while (<FILE>) { 
    push(@linesInFile, $_);
}

close FILE;


my @fixedproddec = FixProductionsAndDecays(\@linesInFile);
my @fixedisos = FixIsos(\@linesInFile);
my @fixedphase = Fixpbarp(\@linesInFile);
my @fixedpipiswave = FixPiPiSWave(\@linesInFile);
my @fixedkpiswave = FixKPiSWaveI12(\@linesInFile);
push(@fixedkpiswave, FixKPiSWaveI32(\@linesInFile));
my @fixedgenerickmatrix = FixGenericKMatrix(\@linesInFile);

my $numfixedproddec = @fixedproddec;
my $numfixedphase = @fixedphase;
my $numfixedpipiswave = @fixedpipiswave;
my $numfixedkpiswave = @fixedkpiswave;
my $numfixedisos = @fixedisos;
my $numfixedgenerickmatrix = @fixedgenerickmatrix;

my $numfixed = $numfixedproddec + $numfixedphase + $numfixedgenerickmatrix + 
               $numfixedpipiswave + $numfixedkpiswave + $numfixedisos;

print "\n\n";
print "##################################################\n";
print "##                                              ##\n";
print "##  Output from the pawianfixparameters script  ##\n";
print "##                                              ##\n";
print "##################################################\n";
print "#\n"; 
print "# Fixing $numfixed parameters\n";
print "#\n";
print "# !! The following list is given without any        !!\n";
print "# !! guarantee and has to be checked by the analyst !!\n";
print "#\n\n";

my $fixedline;

print "##################################################\n";
print "# Fixing one pbarp singlet and one triplet phase #\n";
print "##################################################\n";
foreach $fixedline (@fixedphase){
    print "mnParFix = ".$fixedline."\n";
}

print "\n";
print "##################################################\n";
print "# Fixing $numfixedproddec production and decay parameters\n";
print "##################################################\n";
foreach $fixedline (@fixedproddec){
    print "mnParFix = ".$fixedline."\n";
}

print "\n";
print "##################################################\n";
print "# Fixing $numfixedisos isospin parameters\n";
print "##################################################\n";
foreach $fixedline (@fixedisos){
    print "mnParFix = ".$fixedline."\n";
}

print "\n";
print "##################################################\n";
print "# Fixing $numfixedpipiswave PiPiS-wave parameters\n";
print "##################################################\n";
foreach $fixedline (@fixedpipiswave){
    print "mnParFix = ".$fixedline."\n";
}

print "\n";
print "##################################################\n";
print "# Fixing $numfixedkpiswave KPiS-wave parameters\n";
print "##################################################\n";
foreach $fixedline (@fixedkpiswave){
    print "mnParFix = ".$fixedline."\n";
}

print "\n";
print "##################################################\n";
print "# Fixing $numfixedgenerickmatrix generic k-matrix parameters\n";
print "##################################################\n";
foreach $fixedline (@fixedgenerickmatrix){
    print "mnParFix = ".$fixedline."\n";
}



sub FixProductionsAndDecays
{
    my @linesInFile = @{$_[0]};
    my @fixedParams;
    my @fixedSources;

    foreach my $currentLine (@linesInFile){

	if($currentLine =~ /Matrix/){
	    next;
	} 

	if($currentLine =~ /(.*_(.*)To(.*))Mag.*/ ){
	    my $currentSource = $2;
	    my $currentDest = $3;
	    my $found = 0;

	    if(grep $_ eq $currentSource, @fixedSources){
		$found = 1;
	    }

	    # Fix the priority system instead of this one
	    elsif(!($currentDest eq $prio)){	    
		foreach my $alternative(@linesInFile){
		    if( $alternative =~ /(.*_(.*)To(.*))Mag.*/ ){
			if(($currentSource eq $2) && ($prio eq $3)){
			    $found = 1;
			    last;
			}
		    }
		}
	    }
	    if($found == 0){
		push(@fixedSources, $currentSource);
		push(@fixedParams, $1."Mag");
		push(@fixedParams, $1."Phi");
	    }
	}
    }

    return @fixedParams;
}



sub FixIsos
{
    my @linesInFile = @{$_[0]};
    my @fixedParams;
    my @fixedSources;

    foreach my $currentLine (@linesInFile){
	if( $currentLine =~ /(Iso.+)\s+(\S+)\s+\S+/){
	    push(@fixedParams, $1);
	}
    }

    return @fixedParams;
}



sub Fixpbarp
{
    my @linesInFile = @{$_[0]};
    my @fixedParams;
    my @fixedSources;

    my $pbarpSingletIsFixed=0;
    my $pbarpTripletIsFixed=0;

    foreach my $currentLine (@linesInFile){

	# Fix pbarp phases
	if($currentLine =~ /(J.+P(.+)C(.+)L.+S.+pbarpPhi).*/){
	    if((($2 eq "1" && $3 eq "1") || ($2 eq "-1" && $3 eq "-1")) && $pbarpTripletIsFixed == 0){
		push(@fixedParams, $1);
		$pbarpTripletIsFixed = 1;
	    }
	    elsif((($2 eq "1" && $3 eq "-1") || ($2 eq "-1" && $3 eq "1")) && $pbarpSingletIsFixed == 0){
		push(@fixedParams, $1);
		$pbarpSingletIsFixed = 1;
	    }
	}

    }

    if(!$pbarpSingletIsFixed || !$pbarpTripletIsFixed){
	die "Failed to fix pbarp phases.";
    }

    return @fixedParams;
}



sub FixPiPiSWave
{
    my @linesInFile = @{$_[0]};
    my @fixedParams;
    my @fixedSources;

    foreach my $currentLine (@linesInFile){
	
	if( $currentLine =~ /(.*_(.*)TopipiS.*)Mag.*/ ){
	    my $currentSource = $2;

	    if(!(grep $_ eq $currentSource, @fixedSources)){
		push(@fixedSources, $currentSource);
		push(@fixedParams, $1."Mag");
		push(@fixedParams, $1."Phi");
	    }

	}
	elsif( $currentLine =~ /(pipiS.*S0_PosNegOther).*/){
	    push(@fixedParams, $1);
	}
    }

    return @fixedParams;
}




sub FixKPiSWaveI12
{
    my @linesInFile = @{$_[0]};
    my @fixedParams;

    foreach my $currentLine (@linesInFile){

	if( $currentLine =~ /(.*_(.*)ToKpiS12_.*)Mag.*/ ){
	    my $currentSource = $2;
	    push(@fixedParams, $1."Mag");
	    push(@fixedParams, $1."Phi");
	}
    }
    
    return @fixedParams;
}



sub FixKPiSWaveI32
{
    my @linesInFile = @{$_[0]};
    my @fixedParams;

    foreach my $currentLine (@linesInFile){

	if($currentLine =~ /(.*_(.*)ToKpiS32_.*)Mag.*/ ){
	    my $currentSource = $2;
	    push(@fixedParams, $1."Mag");
	    push(@fixedParams, $1."Phi");
	}
    }
    
    return @fixedParams;
}



sub FixGenericKMatrix
{
    my @linesInFile = @{$_[0]};
    my @fixedParams;
    my @fixedJPCKMat;

    foreach my $currentLine (@linesInFile){

	if($currentLine =~ /(.*Matrix.*To.*)Mag/){
    	    push(@fixedParams, $1."Mag");
	    push(@fixedParams, $1."Phi");
	}
	elsif($currentLine =~ /(.*_(.*To.*Matrix.*))Mag/){
	    if(grep $_ eq $2, @fixedJPCKMat){
		next;
	    }

	    push(@fixedParams, $1."Mag");
	    push(@fixedParams, $1."Phi");
	    push(@fixedJPCKMat, $2);
	}
    }

    return @fixedParams;
}



sub Help
{

    print STDERR <<INLINE_LITERAL_TEXT;

NAME
    pawianfixparameters - automatically fix parameters for PAWIAN Partial Wave Analyses

SYNOPSIS
    pawianfixparameters -f <file> [-p <priority>]

OPTIONS

    -h | --help     Show this help
    -f | --file     Select the file containing the fit parameters
    -p | --prio     Set a priority production system that becomes fixed
                    rather than another system

EXAMPLES

    pawianfixparameters -f defaultparams.dat

    pawianfixparameters -f myparams.dat -prio "K*K" 
    #Fixes K*K amplitudes instead of K*(1430)K amplitudes


INLINE_LITERAL_TEXT
}
