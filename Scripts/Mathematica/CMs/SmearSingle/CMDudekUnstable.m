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

plotCMDudekStable = ReImPlot[cmDudek[(b + I valImag) (b + I valImag), mRec, MBW], {b, 0., 4.0}, PlotStyle -> {Red,{Red,Dashed}}, PlotLegends->Automatic, AxesLabel->{Sqrt[s] [GeV/(c c)], CM Dudek stable}]; 
Export["CMDudekStable.pdf", plotCMDudekStable];

plotCMDudekviaBBUnstable = ReImPlot[RNIntDudekviaBB[(b + I valImag) (b + I valImag)][mDec1, mDec2, mRec, f2PoleBB, m2PoleBB, (MBW MBW), orbMomBW, scalingFacBBSmearing], {b, 0., 4.0}, PlotStyle -> {Blue,{Blue,Dashed}}, PlotLegends->Automatic, AxesLabel->{Sqrt[s] [GeV/(c c)], CM Dudek via BB unstable}];
Export["CMDudekViaBBUnstable.pdf", plotCMDudekviaBBUnstable];

plCMDudekviaBBComp=Show[plotCMDudekStable, plotCMDudekviaBBUnstable];
Export["CMDudekViaBB_compare.pdf", plCMDudekviaBBComp];

 plotCMDudekviaDudekUnstable = ReImPlot[RNIntDudekviaDudek[(b + I valImag) (b + I valImag)][mDec1, mDec2, mRec, f2PoleDudek, m2PoleDudek, (MBW MBW), orbMomBW, scalingFacDudekSmearing], {b, 0., 4.0}, PlotStyle -> {Purple,{Purple,Dashed}}, PlotLegends->Automatic, AxesLabel->{Sqrt[s] [GeV/(c c)], CM Dudek via Dudek unstable}];
Export["CMDudekViaDudekUnstable.pdf", plotCMDudekviaDudekUnstable];

plCMDudekviaDudekComp=Show[plotCMDudekStable, plotCMDudekviaDudekUnstable];
Export["CMDudekViaDudek_compare.pdf", plCMDudekviaDudekComp];


plCMDudekviaDudekvsBBComp=Show[plotCMDudekviaBBUnstable, plotCMDudekviaDudekUnstable];
Export["CMDudekviaDudekvsBB_compare.pdf", plCMDudekviaDudekvsBBComp];


plotCMDudekviaBWUnstable = ReImPlot[RNIntDudekviaBW[(b + I valImag) (b + I valImag)][mDec1, mDec2, mRec, orbMomBW, MBW, GammaBW, scalingFacBWSmearing], {b, 0., 4.0}, PlotStyle -> {Purple,{Purple,Dashed}}, PlotLegends->Automatic, AxesLabel->{Sqrt[s] [GeV/(c c)], CM Dudek via BW unstable}];
Export["CMDudekViaBWUnstable.pdf", plotCMDudekviaBWUnstable];

plCMDudekviaBWComp=Show[plotCMDudekStable, plotCMDudekviaBWUnstable];
Export["CMDudekViaBW_compare.pdf", plCMDudekviaBWComp];

plCMDudekviaDudekvsBWComp=Show[plotCMDudekviaBWUnstable, plotCMDudekviaDudekUnstable];
Export["CMDudekviaDudekvsBW_compare.pdf", plCMDudekviaDudekvsBWComp];
