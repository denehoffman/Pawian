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
Get[FileNameJoin[{Directory[], "../InputParamsDoubleSigmaSigma.m"}]]

(* Precision[1.2]; *)
(* Precision[12/10]; *)

Print["RNIntDudekDoubleViaDudek[(1.0 - I 0.1) (1.0 - I 0.1)][mDec11, mDec12, f2PoleDudek1, m2PoleDudek1, orbMomBW1, (mDec11+mDec12)^2, mDec21, mDec22, f2PoleDudek2, m2PoleDudek2, orbMomBW2,(mDec21+mDec22)^2 , scalingFacDudekSmearing]: " , RNIntDudekDoubleViaDudek[(1.0 - I 0.1) (1.0 - I 0.1)][mDec11, mDec12, f2PoleDudek1, m2PoleDudek1, orbMomBW1, (mDec11+mDec12)^2, mDec21, mDec22, f2PoleDudek2, m2PoleDudek2, orbMomBW2,(mDec21+mDec22)^2 , scalingFacDudekSmearing]];
      
plotCMDudekDoubleviaDudekUnstable = ReImPlot[RNIntDudekDoubleViaDudek[(b + I valImag) (b + I valImag)][mDec11, mDec12, f2PoleDudek1, m2PoleDudek1, orbMomBW1, (mDec11+mDec12)^2, mDec21, mDec22, f2PoleDudek2, m2PoleDudek2, orbMomBW2,(mDec21+mDec22)^2 , scalingFacDudekSmearing], {b, 0., 2.5}, PlotStyle -> {Purple,{Purple,Dashed}}, PlotLegends->Automatic, AxesLabel->{Sqrt[s] [GeV/(c c)], CM Dudek double via Dudek unstable}];

Export["CMDudekDoubleViaDudek.pdf", plotCMDudekDoubleviaDudekUnstable];

