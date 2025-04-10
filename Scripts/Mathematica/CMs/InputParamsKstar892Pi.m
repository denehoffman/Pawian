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
mK:=0.49761
mpi:=0.13957;
mpi0:=0.1349768;

valImag := 0.00;

(**** parameters K*(892)0 pi0 -> (K- pi+) pi0 based on BB CM ****)
(* Mass of the unstable resonance*)
MBW:=0.8955;

(*Width of unstable resonance*)
GammaBW:=0.0514;

(* Masses of the decay particles *)
mDec1:=mK;
mDec2:=mpi;

(* Mass of the decay recoil particle *)
mRec=mpi0;

orbMomBW = 1;


f2PoleBB := 0.31133809685087025051448922330820490963;
m2PoleBB := 0.8238329330807920326010290727540961481;

f2PoleDudek := 0.28759507372364078839913006596975963494;
m2PoleDudek := 0.78197101092932291878039859817243872515;

scalingFacBBSmearing:= 1. / 1.105361310555928;
scalingFacDudekSmearing:=1. / 1.105361310555928;
scalingFacBWSmearing:= 1.;
(* scalingFacBWSmearing:= 1. / 0.08802808902167991; *)
