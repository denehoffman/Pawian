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

(**** parameters rho(700)+ pi- -> (pi+ pi0) pi- based on BB CM ****)
(* Mass of the unstable resonance*)
MBW:=0.77511;

(*Width of unstable resonance*)
GammaBW:=0.1491;

(* Masses of the decay particles *)
mDec1:=mpi+mpi0;
mDec2:=mpi;

(* Mass of the decay recoil particle *)
mRec=mpi0;

orbMomBW = 1;


f2PoleBB := 0.4207051462764408283327900815775571806;
m2PoleBB := 0.60605788493333204795276935269551359375;

f2PoleDudek := 0.36966524147228358938226916212290303466;
m2PoleDudek := 0.52303101506895045228388878328917860629;

scalingFacBBSmearing:= 1./ 1.0655618168689487;
scalingFacDudekSmearing:=1. / 1.1786373003089863;
scalingFacBWSmearing:= 1. / 0.7697543609311899;
(* scalingFacBWSmearing:= 1. / 0.08802808902167991; *)
