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
(* Get[FileNameJoin[{Directory[], "../InputParamsRho770Pi.m"}]] *)

Get[FileNameJoin[{Directory[], "../InputParamsKstar892Pi.m"}]]
(* Precision[1.2]; *)
(* Precision[12/10]; *)

plotIntegrFuncDudek3Mass = Plot[ {f2PoleDudek ImagSigmads2Dudek[b b, mDec1, mDec2, f2PoleDudek, m2PoleDudek, (MBW MBW), orbMomBW]},{b, (mDec1+mDec2), 3.}, AxesLabel->{M [GeV/(c c)], integr} ];
Export["IntegrFuncDudek3Mass.pdf", plotIntegrFuncDudek3Mass];


plotIntegrFuncDudek3 = Plot[ {f2PoleDudek ImagSigmads2Dudek[b, mDec1, mDec2, f2PoleDudek, m2PoleDudek, (MBW MBW), orbMomBW]},{b, (mDec1+mDec2)^2, 3.}, AxesLabel->{sprime [GeV/(c c)], integr} ];
Export["IntegrFuncDudek3.pdf", plotIntegrFuncDudek3];

plotIntegrFuncDudek10 = Plot[ {f2PoleDudek ImagSigmads2Dudek[b, mDec1, mDec2, f2PoleDudek, m2PoleDudek, (MBW MBW), orbMomBW]},{b, (mDec1+mDec2)^2, 10.}, AxesLabel->{sprime [GeV/(c c c c)], integr} ];
Export["IntegrFuncDudek10.pdf", plotIntegrFuncDudek10];

plotIntegrFuncDudek100 = Plot[ {f2PoleDudek ImagSigmads2Dudek[b, mDec1, mDec2, f2PoleDudek, m2PoleDudek, (MBW MBW), orbMomBW]},{b, (mDec1+mDec2)^2, 100.}, AxesLabel->{sprime [GeV GeV/(c c c c)], integr} ];
Export["IntegrFuncDudek100.pdf", plotIntegrFuncDudek100];

plotIntegrFuncDudek1000 = Plot[ {f2PoleDudek ImagSigmads2Dudek[b, mDec1, mDec2, f2PoleDudek, m2PoleDudek, (MBW MBW), orbMomBW]},{b, (mDec1+mDec2)^2, 1000.}, AxesLabel->{sprime [GeV GeV/(c c c c)], integr} ];
Export["IntegrFuncDudek1000.pdf", plotIntegrFuncDudek1000];

plotIntegrFuncDudek10000 = Plot[ {f2PoleDudek ImagSigmads2Dudek[b, mDec1, mDec2, f2PoleDudek, m2PoleDudek, (MBW MBW), orbMomBW]},{b, (mDec1+mDec2)^2, 10000.}, AxesLabel->{sprime [GeV GeV/(c c c c)], integr} ];
Export["IntegrFuncDudek10000.pdf", plotIntegrFuncDudek10000];

Print["intergal s<3: ", RNIntDudek[mDec1, mDec2, f2PoleDudek, m2PoleDudek, (MBW MBW), orbMomBW, scalingFacDudekSmearing, 3.]];
Print["intergal s<100000: ", RNIntDudek[mDec1, mDec2, f2PoleDudek, m2PoleDudek, (MBW MBW), orbMomBW, scalingFacDudekSmearing, 100000.]];

 plotIntegralDudeksmax10 = Plot[{RNIntDudek[mDec1, mDec2, f2PoleDudek, m2PoleDudek, (MBW MBW), orbMomBW, scalingFacDudekSmearing, b]}, {b, 3., 10.}, AxesLabel->{smax [GeV GeV/(c c c c)], integral} ]; 
Export["IntegralDudeksmax10.pdf", plotIntegralDudeksmax10];

plotIntegralDudeksmax100 = Plot[{RNIntDudek[mDec1, mDec2, f2PoleDudek, m2PoleDudek, (MBW MBW), orbMomBW, scalingFacDudekSmearing, b]}, {b, 3., 100.}, AxesLabel->{smax [GeV GeV/(c c c c)], integral} ]; 
Export["IntegralDudeksmax100.pdf", plotIntegralDudeksmax100];

plotIntegralDudeksmax1000 = Plot[{RNIntDudek[mDec1, mDec2, f2PoleDudek, m2PoleDudek, (MBW MBW), orbMomBW, scalingFacDudekSmearing, b]}, {b, 3., 1000.}, AxesLabel->{smax [GeV GeV/(c c c c)], integral} ]; 
Export["IntegralDudeksmax1000.pdf", plotIntegralDudeksmax1000];

plotIntegralDudeksmax10000 = Plot[{RNIntDudek[mDec1, mDec2, f2PoleDudek, m2PoleDudek, (MBW MBW), orbMomBW, scalingFacDudekSmearing, b]}, {b, 3., 10000.}, AxesLabel->{smax [GeV GeV/(c c c c)], integral} ]; 
Export["IntegralDudeksmax10000.pdf", plotIntegralDudeksmax10000];
 
