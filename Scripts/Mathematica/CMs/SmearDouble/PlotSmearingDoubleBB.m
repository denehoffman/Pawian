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

plotIntegrFuncBB3Mass1 = Plot[ {f2PoleBB1 ImagSigmads2BB[b b, mDec11, mDec12, f2PoleBB1, m2PoleBB1, (MBW1 MBW1), orbMomBW1]},{b, (mDec11+mDec12), 3.}, AxesLabel->{M [GeV/(c c)], integr} ];
Export["IntegrFuncBB3Mass1.pdf", plotIntegrFuncBB3Mass1];

plotIntegrFuncBB3Mass2 = Plot[ {f2PoleBB2 ImagSigmads2BB[b b, mDec21, mDec22, f2PoleBB2, m2PoleBB2, (MBW2 MBW2), orbMomBW2]},{b, (mDec21+mDec22), 3.}, AxesLabel->{M [GeV/(c c)], integr} ];
Export["IntegrFuncBB3Mass2.pdf", plotIntegrFuncBB3Mass2];

Print["RNIntBBDouble[mDec11, mDec12, f2PoleBB1, m2PoleBB1, orbMomBW1, (mDec11+mDec12)^2, mDec21, mDec22, f2PoleBB2, m2PoleBB2, orbMomBW2,(mDec21+mDec22)^2 , scalingFacBBSmearing, 10.]: ", RNIntBBDouble[mDec11, mDec12, f2PoleBB1, m2PoleBB1, orbMomBW1, (mDec11+mDec12)^2, mDec21, mDec22, f2PoleBB2, m2PoleBB2, orbMomBW2,(mDec21+mDec22)^2 , scalingFacBBSmearing, 10.]];

Print["RNIntBBDouble[mDec11, mDec12, f2PoleBB1, m2PoleBB1, orbMomBW1, (mDec11+mDec12)^2, mDec21, mDec22, f2PoleBB2, m2PoleBB2, orbMomBW2,(mDec21+mDec22)^2 , scalingFacBBSmearing, 1000.]: ", RNIntBBDouble[mDec11, mDec12, f2PoleBB1, m2PoleBB1, orbMomBW1, (mDec11+mDec12)^2, mDec21, mDec22, f2PoleBB2, m2PoleBB2, orbMomBW2,(mDec21+mDec22)^2 , scalingFacBBSmearing, 1000.]];

Print["RNIntBBDouble[mDec11, mDec12, f2PoleBB1, m2PoleBB1, orbMomBW1, (mDec11+mDec12)^2, mDec21, mDec22, f2PoleBB2, m2PoleBB2, orbMomBW2,(mDec21+mDec22)^2 , scalingFacBBSmearing, 100000.]: ", RNIntBBDouble[mDec11, mDec12, f2PoleBB1, m2PoleBB1, orbMomBW1, (mDec11+mDec12)^2, mDec21, mDec22, f2PoleBB2, m2PoleBB2, orbMomBW2,(mDec21+mDec22)^2 , scalingFacBBSmearing, 100000.]];


