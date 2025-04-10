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

(* Precision[12/10]; *)

(* masses of the decay products *)
mpi:=0.13957;
mpi0:=0.1349768;

mKp:=0.493677;
mKm:=0.493677;

valImag := 0.00;

(* Masses of the unstable resonance*)
MBW1:=0.78266;
MBW2:=1.019461;

(*Widths of unstable resonance*)
GammaBW1:=0.00868;
GammaBW2:=0.004249;

(* Masses of the decay particles *)
mDec11:=2 mpi;
mDec12:=mpi0;

mDec21:=mKp;
mDec22:=mKm;


orbMomBW1 = 1;
orbMomBW2 = 1;

f2PoleBB1 := 0.02275031992477571786636110849745549994;
m2PoleBB1 := 0.61238097137475222594286512618182349074;

f2PoleBB2 := 0.67110752879828435237795796451561483007;
m2PoleBB2 := 1.05187282169299207701471933416642866246;

f2PoleDudek1 := 0.0225865133682695854973481933446531801;
m2PoleDudek1 := 0.60741971755426057356854570645121124506;

f2PoleDudek2 := 0.5596430172817635954719483569269226501;
m2PoleDudek2 := 1.03873920225401985970974681359698101954;


scalingFacBBSmearing:= 1. / 3.287318346901323;
scalingFacDudekSmearing:=1. / 3.8866039897112796;
scalingFacBWSmearing:= 1. / 0.7454720957505188;

