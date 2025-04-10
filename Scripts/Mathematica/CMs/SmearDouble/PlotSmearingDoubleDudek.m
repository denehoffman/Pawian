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

plotIntegrFuncDudek3Mass1 = Plot[ {f2PoleDudek1 ImagSigmads2Dudek[b b, mDec11, mDec12, f2PoleDudek1, m2PoleDudek1, (MBW1 MBW1), orbMomBW1]},{b, (mDec11+mDec12), 3.}, AxesLabel->{M [GeV/(c c)], integr} ];
Export["IntegrFuncDudek3Mass1.pdf", plotIntegrFuncDudek3Mass1];

plotIntegrFuncDudek3Mass2 = Plot[ {f2PoleDudek2 ImagSigmads2Dudek[b b, mDec21, mDec22, f2PoleDudek2, m2PoleDudek2, (MBW2 MBW2), orbMomBW2]},{b, (mDec21+mDec22), 3.}, AxesLabel->{M [GeV/(c c)], integr} ];
Export["IntegrFuncDudek3Mass2.pdf", plotIntegrFuncDudek3Mass2];

Print["RNIntDudekDouble[mDec11, mDec12, f2PoleDudek1, m2PoleDudek1, orbMomBW1, (mDec11+mDec12)^2, mDec21, mDec22, f2PoleDudek2, m2PoleDudek2, orbMomBW2,(mDec21+mDec22)^2 , scalingFacDudekSmearing, 10.]: ", RNIntDudekDouble[mDec11, mDec12, f2PoleDudek1, m2PoleDudek1, orbMomBW1, (mDec11+mDec12)^2, mDec21, mDec22, f2PoleDudek2, m2PoleDudek2, orbMomBW2,(mDec21+mDec22)^2 , scalingFacDudekSmearing, 10.]];

Print["RNIntDudekDouble[mDec11, mDec12, f2PoleDudek1, m2PoleDudek1, orbMomBW1, (mDec11+mDec12)^2, mDec21, mDec22, f2PoleDudek2, m2PoleDudek2, orbMomBW2,(mDec21+mDec22)^2 , scalingFacDudekSmearing, 1000.]: ", RNIntDudekDouble[mDec11, mDec12, f2PoleDudek1, m2PoleDudek1, orbMomBW1, (mDec11+mDec12)^2, mDec21, mDec22, f2PoleDudek2, m2PoleDudek2, orbMomBW2,(mDec21+mDec22)^2 , scalingFacDudekSmearing, 1000.]];

Print["RNIntDudekDouble[mDec11, mDec12, f2PoleDudek1, m2PoleDudek1, orbMomBW1, (mDec11+mDec12)^2, mDec21, mDec22, f2PoleDudek2, m2PoleDudek2, orbMomBW2,(mDec21+mDec22)^2 , scalingFacDudekSmearing, 100000.]: ", RNIntDudekDouble[mDec11, mDec12, f2PoleDudek1, m2PoleDudek1, orbMomBW1, (mDec11+mDec12)^2, mDec21, mDec22, f2PoleDudek2, m2PoleDudek2, orbMomBW2,(mDec21+mDec22)^2 , scalingFacDudekSmearing, 100000.]];


