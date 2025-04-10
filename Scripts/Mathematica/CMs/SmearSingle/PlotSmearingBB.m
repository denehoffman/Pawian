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


plotIntegrFuncBB3Mass = Plot[ {f2PoleBB ImagSigmads2BB[b b, mDec1, mDec2, f2PoleBB, m2PoleBB, (MBW MBW), orbMomBW]},{b, (mDec1+mDec2), 3.}, AxesLabel->{M [GeV/(c c)], integr} ];
Export["IntegrFuncBB3Mass.pdf", plotIntegrFuncBB3Mass];

plotIntegrFuncBB3 = Plot[ {f2PoleBB ImagSigmads2BB[b, mDec1, mDec2, f2PoleBB, m2PoleBB, (MBW MBW), orbMomBW]},{b, (mDec1+mDec2)^2, 3.}, AxesLabel->{sprime [GeV/(c c)], integr} ];
Export["IntegrFuncBB3.pdf", plotIntegrFuncBB3];

plotIntegrFuncBB10 = Plot[ {f2PoleBB ImagSigmads2BB[b, mDec1, mDec2, f2PoleBB, m2PoleBB, (MBW MBW), orbMomBW]},{b, (mDec1+mDec2)^2, 10.}, AxesLabel->{sprime [GeV/(c c c c)], integr} ];
Export["IntegrFuncBB10.pdf", plotIntegrFuncBB10];

plotIntegrFuncBB100 = Plot[ {f2PoleBB ImagSigmads2BB[b, mDec1, mDec2, f2PoleBB, m2PoleBB, (MBW MBW), orbMomBW]},{b, (mDec1+mDec2)^2, 100.}, AxesLabel->{sprime [GeV GeV/(c c c c)], integr} ];
Export["IntegrFuncBB100.pdf", plotIntegrFuncBB100];

plotIntegrFuncBB1000 = Plot[ {f2PoleBB ImagSigmads2BB[b, mDec1, mDec2, f2PoleBB, m2PoleBB, (MBW MBW), orbMomBW]},{b, (mDec1+mDec2)^2, 1000.}, AxesLabel->{sprime [GeV GeV/(c c c c)], integr} ];
Export["IntegrFuncBB1000.pdf", plotIntegrFuncBB1000];

plotIntegrFuncBB10000 = Plot[ {f2PoleBB ImagSigmads2BB[b, mDec1, mDec2, f2PoleBB, m2PoleBB, (MBW MBW), orbMomBW]},{b, (mDec1+mDec2)^2, 10000.}, AxesLabel->{sprime [GeV GeV/(c c c c)], integr} ];
Export["IntegrFuncBB10000.pdf", plotIntegrFuncBB10000];

Print["intergal s<3: ", RNIntBB[mDec1, mDec2, f2PoleBB, m2PoleBB, (MBW MBW), orbMomBW, scalingFacBBSmearing, 3.]];
Print["intergal s<100000: ", RNIntBB[mDec1, mDec2, f2PoleBB, m2PoleBB, (MBW MBW), orbMomBW, scalingFacBBSmearing, 100000.]];

plotIntegralBBsmax10 = Plot[{RNIntBB[mDec1, mDec2, f2PoleBB, m2PoleBB, (MBW MBW), orbMomBW, scalingFacBBSmearing, b]}, {b, 3., 10.}, AxesLabel->{smax [GeV GeV/(c c c c)], integral} ]; 
Export["IntegralBBsmax10.pdf", plotIntegralBBsmax10];

plotIntegralBBsmax100 = Plot[{RNIntBB[mDec1, mDec2, f2PoleBB, m2PoleBB, (MBW MBW), orbMomBW, scalingFacBBSmearing, b]}, {b, 3., 100.}, AxesLabel->{smax [GeV GeV/(c c c c)], integral} ]; 
Export["IntegralBBsmax100.pdf", plotIntegralBBsmax100];

plotIntegralBBsmax1000 = Plot[{RNIntBB[mDec1, mDec2, f2PoleBB, m2PoleBB, (MBW MBW), orbMomBW, scalingFacBBSmearing, b]}, {b, 3., 1000.}, AxesLabel->{smax [GeV GeV/(c c c c)], integral} ]; 
Export["IntegralBBsmax1000.pdf", plotIntegralBBsmax1000];

plotIntegralBBsmax10000 = Plot[{RNIntBB[mDec1, mDec2, f2PoleBB, m2PoleBB, (MBW MBW), orbMomBW, scalingFacBBSmearing, b]}, {b, 3., 10000.}, AxesLabel->{smax [GeV GeV/(c c c c)], integral} ]; 
Export["IntegralBBsmax10000.pdf", plotIntegralBBsmax10000];
 
