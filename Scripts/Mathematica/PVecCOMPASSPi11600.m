(*  Copyright 2023                                                                                
 Bertram Kopf (bertram@ep1.rub.de)
 Meike Kuessner (mkuessner@ep1.rub.de)
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
 This Mathematica-script calculates some relevant quantities for the production process (pi p (COMPASS) )based on the K-matrix parametrization of the pi1-wave published in Eur. Phys. J. C81 (2021) no.12, 1056; doi:10.1140/epjc/s10052-021-09821-2; arXiv:2008.11566.
 *)

Get[FileNameJoin[{Directory[], "KmatUtils.m"}]]

(* Precision[1.2]; *)
(* Precision[12/10]; *)

(* masses of the decay products *)
mpi:=0.13957;
meta:=0.547862;
metaprime:=0.95778;

(* K-matrix parameter *)
mpi11600:=1.385518356865105;
gpieta:= 0.8056407725708118;
gpietaprime:= 1.046952875320115;
c00:= 1.0500;
c01:= 0.1516321466620786;
c11:= -0.2461137125335851;

(* P-vector parameter for COMPASS data *)
bCOMPASSpi1Mag:= 101.3643146415384; 
bCOMPASSpi1Phi:= 0.1965729386642642;
bgPVecCOMPASS00pi1:= -7.689428940848806;
bgPVecCOMPASS01pi1:= -139.9999961614917;
bgPVecCOMPASS10pi1:= 14.53171845010211;
bgPVecCOMPASS11pi1:= 39.44684637184774;
pietaprimeCOMPASSScaling := Sqrt[0.4227529865757182];

bCOMPASSpi1Complex := bCOMPASSpi1Mag  ( Cos[bCOMPASSpi1Phi] + I Sin[bCOMPASSpi1Phi] );

q0PiEta:=BreakupMomQDefaultFromS[mpi11600 mpi11600,mpi mpi, meta meta]; 
q0PiEtaprime:=BreakupMomQDefaultFromS[mpi11600 mpi11600,mpi mpi, metaprime metaprime]; 


KmatPi1BW[s_]:=
  Module[{resultMatr},
	 qPiEta = BreakupMomQDefaultFromS[s,mpi mpi, meta meta];
	 qPiEtaprime = BreakupMomQDefaultFromS[s,mpi mpi, metaprime metaprime];
	 BWPiEta = BWBarrier[1, qPiEta, q0PiEta];
	 BWPiEtaprime = BWBarrier[1, qPiEtaprime, q0PiEtaprime]; 
	 resultMatr = {{BWPiEta (gpieta gpieta/(mpi11600 mpi11600-s)+c00) BWPiEta, BWPiEta (gpieta gpietaprime/(mpi11600*mpi11600-s)+c01) BWPiEtaprime},
		       {BWPiEtaprime (gpietaprime gpieta/(mpi11600 mpi11600-s)+c01) BWPiEta, BWPiEtaprime (gpietaprime gpietaprime/(mpi11600 mpi11600-s)+c11) BWPiEtaprime}};
	 resultMatr
	 ];
         
PVecCOMPASSPi1[s_]:=
  Module[{resultMatr},
	 qPiEta = BreakupMomQDefaultFromS[s,mpi mpi, meta meta];
	 qPiEtaprime = BreakupMomQDefaultFromS[s,mpi mpi, metaprime metaprime];
	 BWPiEta = BWBarrier[1, qPiEta, q0PiEta];
	 BWPiEtaprime = BWBarrier[1, qPiEtaprime, q0PiEtaprime]; 
	 resultVec = {{ (bCOMPASSpi1Complex gpieta/(mpi11600 mpi11600-s)+bgPVecCOMPASS00pi1+bgPVecCOMPASS10pi1 s) BWPiEta},
		       {(bCOMPASSpi1Complex gpietaprime/(mpi11600*mpi11600-s)+bgPVecCOMPASS01pi1+bgPVecCOMPASS11pi1 s) BWPiEtaprime}};
	 resultVec
	 ];

(* Print["PVecCOMPASSPi1[1]: ", PVecCOMPASSPi1[1.] // MatrixForm ]; *)

ChewMmat[s_]:={{c[s,mpi,meta], 0.},{0., c[s,mpi,metaprime]}};

KChewMmat[s_]:=KmatPi1BW[s].ChewMmat[s];

Idmat:=IdentityMatrix[2];      
IKChewMmat[s_]:= Idmat + KChewMmat[s];

invIKChewMmat[s_]:= Inverse[IKChewMmat[s]];

Fvec[s_]:=invIKChewMmat[s].PVecCOMPASSPi1[s];
FPiEtaSqr=Plot[ BreakupMomQDefaultFromS[m m,mpi mpi, meta meta] Norm[Fvec[m m + I 0.000001][[1]]] Norm[Fvec[m m + I 0.000001][[1]]], {m, 0.75, 2.0}, AxesLabel->{M [GeV/(c c)], Abs[F(\[Pi]\[Eta])] Abs[F(\[Pi]\[Eta])]}]; 
Export["FPiEtaSqr.pdf", FPiEtaSqr];

FPiEtaprimeSqr=Plot[ pietaprimeCOMPASSScaling BreakupMomQDefaultFromS[m m,mpi mpi, metaprime metaprime] Norm[Fvec[m m + I 0.000001][[2]]] Norm[Fvec[m m + I 0.000001][[2]]], {m, mpi+metaprime, 2.0}, AxesLabel->{M [GeV/(c c)], Abs[F(\[Pi]\[Eta]')] Abs[F(\[Pi]\[Eta]')]}];
Export["FPiEtaprimeSqr.pdf", FPiEtaprimeSqr];

