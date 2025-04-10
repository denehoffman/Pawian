(*  Copyright 2025                                                                                
 Bertram Kopf (bertram@ep1.rub.de)
    Ruhr-Universit<C3><A4>t Bochum

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

***********************************************************************************************
 *)
Get[FileNameJoin[{Directory[], "../KinUtils.m"}]]
Get[FileNameJoin[{Directory[], "../Integrals.m"}]]
(* Get[FileNameJoin[{Directory[], "../InputParamsPi1300Pi.m"}]] *) 
Get[FileNameJoin[{Directory[], "../InputParamsRho770Pi.m"}]] 

(* Precision[1.2]; *)
(* Precision[12/10]; *)

plotCMBBStable = ReImPlot[cmReid[(b + I valImag) (b + I valImag), mRec, MBW], {b, 0., 4.0}, PlotStyle -> {Red,{Red,Dashed}}, PlotLegends->Automatic, AxesLabel->{Sqrt[s] [GeV/(c c)], CM BB stable}]; 
Export["CMBBStable.pdf", plotCMBBStable];

plotCMBBviaBBUnstable = ReImPlot[RNIntBBviaBB[(b + I valImag) (b + I valImag)][mDec1, mDec2, mRec, f2PoleBB, m2PoleBB, (MBW MBW), orbMomBW, scalingFacBBSmearing], {b, 0., 4.0}, PlotStyle -> {Blue,{Blue,Dashed}}, PlotLegends->Automatic, AxesLabel->{Sqrt[s] [GeV/(c c)], CM BB via BB unstable}];
Export["CMBBViaBBUnstable.pdf", plotCMBBviaBBUnstable];

plCMBBviaBBComp=Show[plotCMBBStable, plotCMBBviaBBUnstable];
Export["CMBBViaBB_compare.pdf", plCMBBviaBBComp];

 plotCMBBviaDudekUnstable = ReImPlot[RNIntBBviaDudek[(b + I valImag) (b + I valImag)][mDec1, mDec2, mRec, f2PoleDudek, m2PoleDudek, (MBW MBW), orbMomBW, scalingFacDudekSmearing], {b, 0., 4.0}, PlotStyle -> {Purple,{Purple,Dashed}}, PlotLegends->Automatic, AxesLabel->{Sqrt[s] [GeV/(c c)], CM BB via BB unstable}];
Export["CMBBViaDudekUnstable.pdf", plotCMBBviaDudekUnstable];

plCMBBviaDudekComp=Show[plotCMBBStable, plotCMBBviaDudekUnstable];
Export["CMBBViaDudek_compare.pdf", plCMBBviaDudekComp];

plCMBBviaDudekvsBBComp=Show[plotCMBBviaBBUnstable, plotCMBBviaDudekUnstable];
Export["CMBBviaDudekvsBB_compare.pdf", plCMBBviaDudekvsBBComp];


plotCMBBviaBWUnstable = ReImPlot[RNIntBBviaBW[(b + I valImag) (b + I valImag)][mDec1, mDec2, mRec, orbMomBW, MBW, GammaBW, scalingFacBWSmearing], {b, 0., 4.0}, PlotStyle -> {Purple,{Purple,Dashed}}, PlotLegends->Automatic, AxesLabel->{Sqrt[s] [GeV/(c c)], CM BB via BB unstable}];
Export["CMBBViaBWUnstable.pdf", plotCMBBviaBWUnstable];

plCMBBviaBWComp=Show[plotCMBBStable, plotCMBBviaBWUnstable];
Export["CMBBViaBW_compare.pdf", plCMBBviaBWComp];

plCMBBviaDudekvsBWComp=Show[plotCMBBviaBWUnstable, plotCMBBviaDudekUnstable];
Export["CMBBviaDudekvsBW_compare.pdf", plCMBBviaDudekvsBWComp];
