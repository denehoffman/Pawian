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

valImag := 0.00;

(* Masses of the unstable resonance*)
MBW1:=0.475;
MBW2:=0.475;

(*Widths of unstable resonance*)
GammaBW1:=0.550;
GammaBW2:=0.550;

(* Masses of the decay particles *)
mDec11:=mpi;
mDec12:=mpi;

mDec21:=mpi0;
mDec22:=mpi0;


orbMomBW1 = 0;
orbMomBW2 = 0;

f2PoleBB1 := 0.53525549369732492365844555006462551715;
m2PoleBB1 := 0.02072166534172146370023379334682632026;

f2PoleBB2 := 0.52537598909653228136099257201590318318;
m2PoleBB2 := 0.01303711145733995864252490339790527583;

f2PoleDudek1 := 0.53525549368142545098154277118417779544;
m2PoleDudek1 := -0.32003256520041713522826607205769323751;

f2PoleDudek2 := 0.52537598910691176848144924060406008584;
m2PoleDudek2 := -0.32142763113897706864276934163696308259;


scalingFacBBSmearing:= 1. / 1.1927781712395298;
scalingFacDudekSmearing:=1. / 1.1927781669814008;
scalingFacBWSmearing:= 1./ 1.1245037531410065;

