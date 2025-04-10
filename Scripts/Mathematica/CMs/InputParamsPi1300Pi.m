(*  Copyright 2023                                                                                
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



(**** CM parameters for pi(1300)0 pi0 w/ pi(1300)0 -> (pi+ pi0) pi- ****)

(* Mass of the unstable resonance*)
MBW:=1.300;

(*Width of unstable resonance*)
GammaBW:=0.400;

(* Masses of the decay particles *)
mDec1:=mpi+mpi0;
mDec2:=mpi;

(* Mass of the decay recoil particle *)
mRec=mpi0;

(* Momentum of the decay pi(1300) -> (pi+ pi0) pi- *)
orbMomBW = 0;

f2PoleBB := 0.62562280402090093702928976311097271904;
m2PoleBB := 1.35581006495341548877557531562277255115;

f2PoleDudek := 0.62562280398394301015550122301071064861;
m2PoleDudek := 0.97194445334541372054996097989249217186;

scalingFacBBSmearing:= 1.;
scalingFacDudekSmearing:= 1.;
scalingFacBWSmearing:= 1. / 0.8336427551137701;
(* scalingFacBWSmearing:= 1. / 0.43630464471242414; *)
