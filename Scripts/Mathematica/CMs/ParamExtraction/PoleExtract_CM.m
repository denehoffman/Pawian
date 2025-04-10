(*  Copyright 2020 
    Meike Kuessner
    mail: mkuessner@ep1.rub.de
    Ruhr-Universität Bochum

    modified by Bertram Kopf (bertram@ep1.rub.de)
 
    This file is part of Pawian.

    Pawian is free software: you can redistribute it and/or modify    
    it under the terms of the GNU General Public License as published by 
    the Free Software Foundation, either version 3 of the License, or     
    (at your option) any later version.                       
                                     
    Pawian is distributed in the hope that it will be useful,         
    but WITHOUT ANY WARRANTY; without even the implied warranty of    
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the     
    GNU General Public License for more details.                  
                                     
    You should have received a copy of the GNU General Public License     
    along with Pawian.  If not, see <http://www.gnu.org/licenses/>

***********************************************************************
This Script computes the bare masses and bare coupling strengths of a resonance
by its physical mass, width and masses of the decay products in case of a
two-body decay. 
The calculation is based on K-matrix descriptions with Chew-Mandelstam functions
 as defined in  Phys.Rev.D 19 (1979) 239 ("BB") and in
 Dudek et al  Phys. Rev. D91, 054008 (2015), Appendix B ("Dudek")
***********************************************************************
*)

Get[FileNameJoin[{Directory[], "../KinUtils.m"}]]
(* Get[FileNameJoin[{Directory[], "../InputParamsPi1300Pi.m"}]] *)
(* Get[FileNameJoin[{Directory[], "../InputParamsRho770Pi.m"}]] *)

Get[FileNameJoin[{Directory[], "../InputParamsKstar892Pi.m"}]]

Print["MBW: ", MBW,"\tGammaBW: ", GammaBW,"\tmDec1: ", mDec1,"\tmDec2: ", mDec2,"\torbL: ", orbMomBW];

(* Singularity of denominator |d(s)|^2 of dispersion relation defines pole position, see Basdevant & Berger, Phys.Rev.D 19 (1979) 239  eq(12) and eq(14) *)

ds2BB[s_,m1_,m2_,fR2_,mR2_, s0_, orbL_]:= Abs[mR2 - s - fR2  BarrierSimple[s, (m1+m2)^2, s0, orbL] Conjugate[cmReid[s,m1,m2]]]^2;
	 
        

(* Determine Minimum and save *)
resultBB=FindMinimum[SetPrecision[ds2BB[(MBW - I GammaBW /2.)^2,mDec1,mDec2,a,b, (MBW MBW), orbMomBW],30],{a,b},WorkingPrecision->20];
{fsqBB,msqBB}={a,b}/.resultBB[[2]]; 
Print["m^2(BB) = ", msqBB, " ", "f^2(BB) = ", fsqBB];
Print["m(BB) = ", Sqrt[msqBB], " ", "f(BB) = ", Sqrt[fsqBB]];

(* Plot Contour around Minimum and Save Plot, change plot style here if needed or comment to not save the plot*)
plotBB=ContourPlot[Log[ds2BB[(MBW-I 2 GammaBW)^2,mDec1,mDec2,f,m, (MBW MBW), orbMomBW]],{f,fsqBB-0.1,fsqBB+0.1},{m,msqBB-0.1,msqBB+0.1},Contours->100,PlotLegends->Automatic,WorkingPrecision->30,ContourStyle->None, ColorFunction->"Rainbow", Epilog->({Black,PointSize[.03],Point[{fsqBB,msqBB}],Text["minimum BB",{fsqBB,msqBB*0.975}]}),FrameLabel->{f^2PoleBB,m^2PoleBB}];

Export["PoleBB.pdf",plotBB];
  

(* Extraction for Dudek CM *)
ds2Dudek[s_,m1_,m2_,fR2_,mR2_, s0_, orbL_]:= Abs[mR2 - s + fR2  BarrierSimple[s, (m1+m2)^2, s0, orbL] Conjugate[cmDudek[s,m1,m2]]]^2;

(* Determine Minimum and save *)
resultDudek=FindMinimum[SetPrecision[ds2Dudek[(MBW - I GammaBW /2.)^2,mDec1,mDec2,a,b, (MBW MBW), orbMomBW],30],{a,b},WorkingPrecision->20];
{fsqDudek,msqDudek}={a,b}/.resultDudek[[2]]; 
Print["m^2(Dudek)= ", msqDudek, " ", "f^2(Dudek)= ", fsqDudek];
Print["m(Dudek)= ", Sqrt[msqDudek], " ", "f(Dudek)= ", Sqrt[fsqDudek]];

(* Plot Contour around Minimum and Save Plot, change plot style here if needed or comment to not save the plot*)
plotDudek=ContourPlot[Log[ds2[(MBW-I 2 GammaBW)^2,mDec1,mDec2,f,m, (MBW MBW), orbMomBW]],{f,fsqDudek-0.1,fsqDudek+0.1},{m,msqDudek-0.1,msqDudek+0.1},Contours->100,PlotLegends->Automatic,WorkingPrecision->30,ContourStyle->None, ColorFunction->"Rainbow", Epilog->({Black,PointSize[.03],Point[{fsqDudek,msqDudek}],Text["minimum Dudek",{fsqDudek,msqDudek*0.975}]}),FrameLabel->{f^2,m^2}];

Export["PoleDudek.pdf",plotDudek];
