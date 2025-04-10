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

plotIntegrFuncBW3Mass = Plot[ {scalingFacBWSmearing rho[b b, mDec1, mDec2] Abs[BreitWignerBlattWRelNorm[b b, orbMomBW, MBW, GammaBW, mDec1, mDec2]]^2},{b, (mDec1+mDec2), 3.}, AxesLabel->{M [GeV/(c c)], integr} ];
Export["IntegrFuncBW3Mass.pdf", plotIntegrFuncBW3Mass];


Print["intergal s<3: ", RNIntBW[MBW, GammaBW, mDec1, mDec2, orbMomBW, scalingFacBWSmearing, 3.]];
Print["intergal s<100000: ", RNIntBW[MBW, GammaBW, mDec1, mDec2, orbMomBW, scalingFacBWSmearing, 100000.]];


plotIntegralBWsmax10 = Plot[{RNIntBW[MBW, GammaBW, mDec1, mDec2, orbMomBW, scalingFacBWSmearing, b]}, {b, 3., 10.}, AxesLabel->{smax [GeV GeV/(c c c c)], integral} ]; 
Export["IntegralBWsmax10.pdf", plotIntegralBWsmax10];

plotIntegralBWsmax100 = Plot[{RNIntBW[MBW, GammaBW, mDec1, mDec2, orbMomBW, scalingFacBWSmearing,b]}, {b, 3., 100.}, AxesLabel->{smax [GeV GeV/(c c c c)], integral} ]; 
Export["IntegralBWsmax100.pdf", plotIntegralBWsmax100];

plotIntegralBWsmax1000 = Plot[{RNIntBW[MBW, GammaBW, mDec1, mDec2, orbMomBW, scalingFacBWSmearing, b]}, {b, 3., 1000.}, AxesLabel->{smax [GeV GeV/(c c c c)], integral} ]; 
Export["IntegralBWsmax1000.pdf", plotIntegralBWsmax1000];

plotIntegralBWsmax10000 = Plot[{RNIntBW[MBW, GammaBW, mDec1, mDec2, orbMomBW, scalingFacBWSmearing, b]}, {b, 3., 10000.}, AxesLabel->{smax [GeV GeV/(c c c c)], integral} ]; 
Export["IntegralBWsmax10000.pdf", plotIntegralBWsmax10000];
 
      
