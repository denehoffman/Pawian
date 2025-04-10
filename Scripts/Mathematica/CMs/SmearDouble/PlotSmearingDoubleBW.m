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
(* Get[FileNameJoin[{Directory[], "../InputParamsDoubleSigmaSigma.m"}]] *)
Get[FileNameJoin[{Directory[], "../InputParamsDoubleOmegaPhi.m"}]]

(* Precision[1.2]; *)
(* Precision[12/10]; *)

plotIntegrFuncBW3Mass1 = Plot[ {rho[b b, mDec11, mDec12] Abs[BreitWignerBlattWRelNorm[b b, orbMomBW1, MBW1, GammaBW1, mDec11, mDec12]]^2},{b, (mDec11+mDec12), 3.}, AxesLabel->{M [GeV/(c c)], integr} ];
Export["IntegrFuncBW3Mass1.pdf", plotIntegrFuncBW3Mass1];

plotIntegrFuncBW3Mass2 = Plot[ {rho[b b, mDec21, mDec22] Abs[BreitWignerBlattWRelNorm[b b, orbMomBW2, MBW2, GammaBW2, mDec21, mDec22]]^2},{b, (mDec21+mDec22), 3.}, AxesLabel->{M [GeV/(c c)], integr} ];
Export["IntegrFuncBW3Mass2.pdf", plotIntegrFuncBW3Mass2];


Print["RNIntBWDouble[MBW1, GammaBW1, mDec11, mDec12, orbMomBW1, MBW2, GammaBW2, mDec21, mDec22, orbMomBW2, scalingFacBWSmearing, 10.]: ", RNIntBWDouble[MBW1, GammaBW1, mDec11, mDec12, orbMomBW1, MBW2, GammaBW2, mDec21, mDec22, orbMomBW2, scalingFacBWSmearing, 10.]];


Print["RNIntBWDouble[MBW1, GammaBW1, mDec11, mDec12, orbMomBW1, MBW2, GammaBW2, mDec21, mDec22, orbMomBW2, scalingFacBWSmearing, 1000.]: ", RNIntBWDouble[MBW1, GammaBW1, mDec11, mDec12, orbMomBW1, MBW2, GammaBW2, mDec21, mDec22, orbMomBW2, scalingFacBWSmearing, 1000.]];

Print["RNIntBWDouble[MBW1, GammaBW1, mDec11, mDec12, orbMomBW1, MBW2, GammaBW2, mDec21, mDec22, orbMomBW2, scalingFacBWSmearing, 100000.]: ", RNIntBWDouble[MBW1, GammaBW1, mDec11, mDec12, orbMomBW1, MBW2, GammaBW2, mDec21, mDec22, orbMomBW2, scalingFacBWSmearing, 100000.]];



