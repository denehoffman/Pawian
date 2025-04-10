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
(* Get[FileNameJoin[{Directory[], "../InputParamsDoubleSigmaSigma.m"}]] *)
Get[FileNameJoin[{Directory[], "../InputParamsDoubleOmegaPhi.m"}]]

Print["MBW1: ", MBW1,"\tGammaBW1: ", GammaBW1,"\tmDec11: ", mDec11,"\tmDec12: ", mDec12,"\torbMomBW1: ", orbMomBW1];
Print["MBW2: ", MBW2,"\tGammaBW2: ", GammaBW2,"\tmDec21: ", mDec21,"\tmDec22: ", mDec22,"\torbMomBW2: ", orbMomBW2];


ds2BB[s_,m1_,m2_,fR2_,mR2_, s0_, orbL_]:= Abs[mR2 - s - fR2  BarrierSimple[s, (m1+m2)^2, s0, orbL] Conjugate[cmReid[s,m1,m2]]]^2;
	 
        

(* Determine Minimum and save *)
resultBB1=FindMinimum[SetPrecision[ds2BB[(MBW1 - I GammaBW1 /2.)^2,mDec11,mDec12,a,b, (MBW1 MBW1), orbMomBW1],30],{a,b},WorkingPrecision->20];
{fsqBB1,msqBB1}={a,b}/.resultBB1[[2]]; 
Print["m^2(BB1) = ", msqBB1, " ", "f^2(BB1) = ", fsqBB1];
Print["m(BB1) = ", Sqrt[msqBB1], " ", "f(BB1) = ", Sqrt[fsqBB1]];

(* Plot Contour around Minimum and Save Plot, change plot style here if needed or comment to not save the plot*)
plotBB1=ContourPlot[Log[ds2BB[(MBW1-I 2 GammaBW1)^2,mDec11,mDec12,f,m, (MBW1 MBW1), orbMomBW1]],{f,fsqBB1-0.1,fsqBB1+0.1},{m,msqBB1-0.1,msqBB1+0.1},Contours->100,PlotLegends->Automatic,WorkingPrecision->30,ContourStyle->None, ColorFunction->"Rainbow", Epilog->({Black,PointSize[.03],Point[{fsqBB,msqBB}],Text["minimum BB 1",{fsqBB1,msqBB1*0.975}]}),FrameLabel->{f^2PoleBB1,m^2PoleBB1}];

Export["PoleBB1.pdf",plotBB1];
  
resultBB2=FindMinimum[SetPrecision[ds2BB[(MBW2 - I GammaBW2 /2.)^2,mDec21,mDec22,a,b, (MBW2 MBW2), orbMomBW2],30],{a,b},WorkingPrecision->20];
{fsqBB2,msqBB2}={a,b}/.resultBB2[[2]]; 
Print["m^2(BB2) = ", msqBB2, " ", "f^2(BB2) = ", fsqBB2];
Print["m(BB2) = ", Sqrt[msqBB2], " ", "f(BB2) = ", Sqrt[fsqBB2]];

(* Plot Contour around Minimum and Save Plot, change plot style here if needed or comment to not save the plot*)
plotBB2=ContourPlot[Log[ds2BB[(MBW2-I 2 GammaBW2)^2,mDec21,mDec22,f,m, (MBW2 MBW2), orbMomBW2]],{f,fsqBB2-0.1,fsqBB2+0.1},{m,msqBB2-0.1,msqBB2+0.1},Contours->100,PlotLegends->Automatic,WorkingPrecision->30,ContourStyle->None, ColorFunction->"Rainbow", Epilog->({Black,PointSize[.03],Point[{fsqBB,msqBB}],Text["minimum BB 2",{fsqBB2,msqBB2*0.975}]}),FrameLabel->{f^2PoleBB2,m^2PoleBB2}];

Export["PoleBB2.pdf",plotBB2];


(* Extraction for Dudek CM *)
ds2Dudek[s_,m1_,m2_,fR2_,mR2_, s0_, orbL_]:= Abs[mR2 - s + fR2  BarrierSimple[s, (m1+m2)^2, s0, orbL] Conjugate[cmDudek[s,m1,m2]]]^2;

(* Determine Minimum and save *)
resultDudek1=FindMinimum[SetPrecision[ds2Dudek[(MBW1 - I GammaBW1 /2.)^2,mDec11,mDec12,a,b, (MBW1 MBW1), orbMomBW1],30],{a,b},WorkingPrecision->20];
{fsqDudek1,msqDudek1}={a,b}/.resultDudek1[[2]]; 
Print["m^2(Dudek1)= ", msqDudek1, " ", "f^2(Dudek1)= ", fsqDudek1];
Print["m(Dudek1)= ", Sqrt[msqDudek1], " ", "f(Dudek1)= ", Sqrt[fsqDudek1]];

(* Plot Contour around Minimum and Save Plot, change plot style here if needed or comment to not save the plot*)
plotDudek1=ContourPlot[Log[ds2[(MBW1-I 2 GammaBW1)^2,mDec11,mDec12,f,m, (MBW1 MBW1), orbMomBW1]],{f,fsqDudek1-0.1,fsqDudek1+0.1},{m,msqDudek1-0.1,msqDudek1+0.1},Contours->100,PlotLegends->Automatic,WorkingPrecision->30,ContourStyle->None, ColorFunction->"Rainbow", Epilog->({Black,PointSize[.03],Point[{fsqDudek1,msqDudek1}],Text["minimum Dudek 1",{fsqDudek1,msqDudek1*0.975}]}),FrameLabel->{f^2,m^2}];

Export["PoleDudek1.pdf",plotDudek1];


resultDudek2=FindMinimum[SetPrecision[ds2Dudek[(MBW2 - I GammaBW2 /2.)^2,mDec21,mDec22,a,b, (MBW2 MBW2), orbMomBW2],30],{a,b},WorkingPrecision->20];
{fsqDudek2,msqDudek2}={a,b}/.resultDudek2[[2]]; 
Print["m^2(Dudek2)= ", msqDudek2, " ", "f^2(Dudek2)= ", fsqDudek2];
Print["m(Dudek2)= ", Sqrt[msqDudek2], " ", "f(Dudek2)= ", Sqrt[fsqDudek2]];

(* Plot Contour around Minimum and Save Plot, change plot style here if needed or comment to not save the plot*)
plotDudek2=ContourPlot[Log[ds2[(MBW2-I 2 GammaBW2)^2,mDec21,mDec22,f,m, (MBW2 MBW2), orbMomBW2]],{f,fsqDudek2-0.1,fsqDudek2+0.1},{m,msqDudek2-0.1,msqDudek2+0.1},Contours->100,PlotLegends->Automatic,WorkingPrecision->30,ContourStyle->None, ColorFunction->"Rainbow", Epilog->({Black,PointSize[.03],Point[{fsqDudek2,msqDudek2}],Text["minimum Dudek 2",{fsqDudek2,msqDudek2*0.975}]}),FrameLabel->{f^2,m^2}];

Export["PoleDudek2.pdf",plotDudek2];
