(*  Copyright 2023                                                                                
 Bertram Kopf (bertram@ep1.rub.de)
 Meike Kuessner (mkuessner@ep1.rub.de)
    Ruhr-Universität Bochum

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
 This Mathematica-script calculates some relevant quantities for the scattering process based on the K-matrix parametrization of the f0-wave published in Eur. Phys. J. C81 (2021) no.12, 1056; doi:10.1140/epjc/s10052-021-09821-2; arXiv:2008.11566.
 *)


Get[FileNameJoin[{Directory[], "KmatUtils.m"}]]

(* Precision[1.2]; *)
  (* Precision[12/10]; *) 

(* masses of the decay products *)
mpi := 0.13957;
(* meta := 0.547; *)
meta := 0.547862;
(* m2pi := 2 mpi; *)
m2pi := 0.26996;
mKp := 0.49367;
mK0 := 0.497614;  
metaprime := 0.95778;

(* K-matrix parameter *)
mf0500 := 0.5146109988244556;
mf0980 := 0.9062999999986513;
mf01370 := 1.23089000002673;
mf01500 := 1.461043944511787;
mf01710 := 1.696114327468766;

gpipif0500 :=  0.749866997688989; 
g4pif0500 := -0.01257099832673861;
gKKf0500 := 0.2753599978535977;
getaetaf0500 := -0.1510199937514032;
getaetaprimef0500 := 0.3610299929020451;

gpipif0980 := 0.06400735441028882; 
g4pif0980 :=  0.002039993700021009;
gKKf0980 := 0.7741299935288173;
getaetaf0980 := 0.5099954460483236;
getaetaprimef0980 := 0.131119996207024;

gpipif01370 :=  -0.2341669602361275;
g4pif01370 := -0.01031664796738707;
gKKf01370 := 0.7228310629513335;
getaetaf01370 := 0.1193373160160431;
getaetaprimef01370 := 0.3679219171982366;

gpipif01500 := 0.01270001206662291;
g4pif01500 := 0.2670000044701449;
gKKf01500 := 0.09214335545338775;
getaetaf01500 := 0.02742288751616556;
getaetaprimef01500 := -0.04024795048926635;

gpipif01710 := -0.1424226773316178;
g4pif01710 := 0.2277971435654336;
gKKf01710 := 0.1598113086438209;
getaetaf01710 := 0.162720778677211;
getaetaprimef01710 := -0.1739657300479793;

c00 := 0.03728069393605827;
c01 := 0.;
c02 := -0.01398000003371962;
c03 := -0.02202999981025169;
c04 := 0.01397000015464572;
c11 := 0.;
c12 := 0.;
c13 := 0.;
c14 := 0.;
c22 := 0.02349000177968514;
c23 := 0.03100999997418123;
c24 := -0.04002991964937379;
c33 := -0.1376928637961125;
c34 := -0.06721849488474475;
c44 := -0.2840099964663654;


s0 := 0.009112500000000001;
s0Adler := 0.1139704455925943;
snormAdler = 1.;
 

Kmatf0[s_]:=
  Module[{resultMatr},
         (* Which[Re[s]<0.01, s=0.1 + I 0.0001]; *)
	 (*	 Print["s: ", s]; *)
 	 
  KPiPitoPiPi:=(gpipif0500 gpipif0500/(mf0500 mf0500-s)+c00)
	 +(gpipif0980 gpipif0980/(mf0980 mf0980-s)+c00)
	 +(gpipif01370 gpipif01370/(mf01370 mf01370-s)+c00)
	 +(gpipif01500 gpipif01500/(mf01500 mf01500-s)+c00)
	 +(gpipif01710 gpipif01710/(mf01710 mf01710-s)+c00);
         (* Print["KPiPitoPiPi: " , KPiPitoPiPi]; *)
	 
  KPiPito4Pi:=(gpipif0500 g4pif0500/(mf0500 mf0500-s)+c01)
         +(gpipif0980 g4pif0980/(mf0980 mf0980-s)+c01)
         +(gpipif01370 g4pif01370/(mf01370 mf01370-s)+c01)
         +(gpipif01500 g4pif01500/(mf01500 mf01500-s)+c01)
         +(gpipif01710 g4pif01710/(mf01710 mf01710-s)+c01);
	 (* Print["KPiPito4Pi: " , KPiPito4Pi];*)
	 
  KPiPitoKK:=(gpipif0500 gKKf0500/(mf0500 mf0500-s)+c02)
         +(gpipif0980 gKKf0980/(mf0980 mf0980-s)+c02)
         +(gpipif01370 gKKf01370/(mf01370 mf01370-s)+c02)
         +(gpipif01500 gKKf01500/(mf01500 mf01500-s)+c02)
         +(gpipif01710 gKKf01710/(mf01710 mf01710-s)+c02);
         (* Print["KPiPitoKK: " , KPiPitoKK]; *)
	 
  KPiPitoEtaEta:=(gpipif0500 getaetaf0500/(mf0500 mf0500-s)+c03)
         +(gpipif0980 getaetaf0980/(mf0980 mf0980-s)+c03)
         +(gpipif01370 getaetaf01370/(mf01370 mf01370-s)+c03)
         +(gpipif01500 getaetaf01500/(mf01500 mf01500-s)+c03)
         +(gpipif01710 getaetaf01710/(mf01710 mf01710-s)+c03);
	 (* Print["KPiPitoEtaEta: " , KPiPitoEtaEta]; *)
	 
         KPiPitoEtaEtaprime:=(gpipif0500 getaetaprimef0500/(mf0500 mf0500-s)+c04)
         +(gpipif0980 getaetaprimef0980/(mf0980 mf0980-s)+c04)
         +(gpipif01370 getaetaprimef01370/(mf01370 mf01370-s)+c04)
         +(gpipif01500 getaetaprimef01500/(mf01500 mf01500-s)+c04)
         +(gpipif01710 getaetaprimef01710/(mf01710 mf01710-s)+c04);
	 (* Print["KPiPitoEtaEtaprime: " , KPiPitoEtaEtaprime]; *)

	 K4Pito4Pi=(g4pif0500 g4pif0500/(mf0500 mf0500-s)+c11)
         +(g4pif0980 g4pif0980/(mf0980 mf0980-s)+c11)
         +(g4pif01370 g4pif01370/(mf01370 mf01370-s)+c11)
         +(g4pif01500 g4pif01500/(mf01500 mf01500-s)+c11)
         +(g4pif01710 g4pif01710/(mf01710 mf01710-s)+c11);

         K4PitoKK=(g4pif0500 gKKf0500/(mf0500 mf0500-s)+c12)
         +(g4pif0980 gKKf0980/(mf0980 mf0980-s)+c12)
         +(g4pif01370 gKKf01370/(mf01370 mf01370-s)+c12)
         +(g4pif01500 gKKf01500/(mf01500 mf01500-s)+c12)
         +(g4pif01710 gKKf01710/(mf01710 mf01710-s)+c12);

         K4PitoEtaEta=(g4pif0500 getaetaf0500/(mf0500 mf0500-s)+c13)
         +(g4pif0980 getaetaf0980/(mf0980 mf0980-s)+c13)
         +(g4pif01370 getaetaf01370/(mf01370 mf01370-s)+c13)
         +(g4pif01500 getaetaf01500/(mf01500 mf01500-s)+c13)
         +(g4pif01710 getaetaf01710/(mf01710 mf01710-s)+c13);
         
         K4PitoEtaEtaprime=(g4pif0500 getaetaprimef0500/(mf0500 mf0500-s)+c14)
         +(g4pif0980 getaetaprimef0980/(mf0980 mf0980-s)+c14)
         +(g4pif01370 getaetaprimef01370/(mf01370 mf01370-s)+c14)
         +(g4pif01500 getaetaprimef01500/(mf01500 mf01500-s)+c14)
         +(g4pif01710 getaetaprimef01710/(mf01710 mf01710-s)+c14);


	 KKKtoKK=(gKKf0500 gKKf0500/(mf0500 mf0500-s)+c22)
         +(gKKf0980 gKKf0980/(mf0980 mf0980-s)+c22)
         +(gKKf01370 gKKf01370/(mf01370 mf01370-s)+c22)
         +(gKKf01500 gKKf01500/(mf01500 mf01500-s)+c22)
         +(gKKf01710 gKKf01710/(mf01710 mf01710-s)+c22);

         KKKtoEtaEta=(gKKf0500 getaetaf0500/(mf0500 mf0500-s)+c23)
         +(gKKf0980 getaetaf0980/(mf0980 mf0980-s)+c23)
         +(gKKf01370 getaetaf01370/(mf01370 mf01370-s)+c23)
         +(gKKf01500 getaetaf01500/(mf01500 mf01500-s)+c23)
         +(gKKf01710 getaetaf01710/(mf01710 mf01710-s)+c23);
         
         KKKtoEtaEtaprime=(gKKf0500 getaetaprimef0500/(mf0500 mf0500-s)+c24)
         +(gKKf0980 getaetaprimef0980/(mf0980 mf0980-s)+c24)
         +(gKKf01370 getaetaprimef01370/(mf01370 mf01370-s)+c24)
         +(gKKf01500 getaetaprimef01500/(mf01500 mf01500-s)+c24)
         +(gKKf01710 getaetaprimef01710/(mf01710 mf01710-s)+c24);


	 KEtaEtatoEtaEta=(getaetaf0500 getaetaf0500/(mf0500 mf0500-s)+c33)
         +(getaetaf0980 getaetaf0980/(mf0980 mf0980-s)+c33)
         +(getaetaf01370 getaetaf01370/(mf01370 mf01370-s)+c33)
         +(getaetaf01500 getaetaf01500/(mf01500 mf01500-s)+c33)
         +(getaetaf01710 getaetaf01710/(mf01710 mf01710-s)+c33);
         
         KEtaEtatoEtaEtaprime=(getaetaf0500 getaetaprimef0500/(mf0500 mf0500-s)+c34)
         +(getaetaf0980 getaetaprimef0980/(mf0980 mf0980-s)+c34)
         +(getaetaf01370 getaetaprimef01370/(mf01370 mf01370-s)+c34)
         +(getaetaf01500 getaetaprimef01500/(mf01500 mf01500-s)+c34)
         +(getaetaf01710 getaetaprimef01710/(mf01710 mf01710-s)+c34);


	 KEtaEtaprimetoEtaEtaprime=(getaetaprimef0500 getaetaprimef0500/(mf0500 mf0500-s)+c44)
         +(getaetaprimef0980 getaetaprimef0980/(mf0980 mf0980-s)+c44)
         +(getaetaprimef01370 getaetaprimef01370/(mf01370 mf01370-s)+c44)
         +(getaetaprimef01500 getaetaprimef01500/(mf01500 mf01500-s)+c44)
         +(getaetaprimef01710 getaetaprimef01710/(mf01710 mf01710-s)+c44);

	adlerTerm=(s - s0)/snormAdler;  
	 

	 resultMatr = {
	   {KPiPitoPiPi,        KPiPito4Pi,        KPiPitoKK,        KPiPitoEtaEta,        KPiPitoEtaEtaprime}, 
	   {KPiPito4Pi,         K4Pito4Pi,         K4PitoKK,         K4PitoEtaEta,         K4PitoEtaEtaprime},
	   {KPiPitoKK,          K4PitoKK,          KKKtoKK,          KKKtoEtaEta,          KKKtoEtaEtaprime},
	   {KPiPitoEtaEta,      K4PitoEtaEta,      KKKtoEtaEta,      KEtaEtatoEtaEta,      KEtaEtatoEtaEtaprime},
	   {KPiPitoEtaEtaprime, K4PitoEtaEtaprime, KKKtoEtaEtaprime, KEtaEtatoEtaEtaprime, KEtaEtaprimetoEtaEtaprime}
	 };
	 resultMatr = adlerTerm resultMatr; 
	 resultMatr
	 ]; 	  
         

ChewMmat[s_]:={{c[s,mpi,mpi], 0., 0., 0., 0.},
	       {0., c[s,m2pi,m2pi], 0. ,0., 0.},
	       {0., 0., c[s,mKp,mK0], 0., 0.},
	       {0., 0., 0., c[s,meta,meta], 0.},
	       {0., 0., 0., 0., c[s,meta,metaprime]}
	 };


KChewMmat[s_]:=Kmatf0[s].ChewMmat[s];
(* Print["CM[1.8]: ", ChewMmat[1.8] // MatrixForm]; *)
      
      
Idmat:=IdentityMatrix[5];      
IKChewMmat[s_]:= Idmat + KChewMmat[s];

invIKChewMmat[s_]:= Inverse[IKChewMmat[s]];
	
Tmat[s_]:=invIKChewMmat[s].Kmatf0[s];

	 (* Print["Tmat[1.8]: ", Tmat[1.8 + I 0.00001] // MatrixForm]; *)    
		  		       
f0PiPitoPiPiTreal=Plot[Re[Tmat[b*b + I 0.000001][[1,1]]],{b, mpi+mpi, 1.9}, AxesLabel->{M [GeV/(c c)], Re[T(\[Pi]\[Pi] -> \[Pi]\[Pi])]}];
Export["f0PiPitoPiPiTreal.pdf", f0PiPitoPiPiTreal];
f0PiPitoPiPiTimag=Plot[Im[Tmat[b*b + I 0.000001][[1,1]]],{b, mpi+mpi, 1.9}, AxesLabel->{M [GeV/(c c)], Im[T(\[Pi]\[Pi] -> \[Pi]\[Pi])]}];
Export["f0PiPitoPiPiTimag.pdf", f0PiPitoPiPiTimag];

f0PiPiElasticity=Plot[Abs[Sii[b*b + I 0.000001, mpi, mpi, Tmat[b*b + I 0.000001][[1,1]]]],{b, mpi+mpi, 1.9}, AxesLabel->{M [GeV/(c c)], \[Eta](\[Pi]\[Pi] -> \[Pi]\[Pi])}];
Export["f0PiPiElasticity.pdf", f0PiPiElasticity];

f0PiPitoPiPiPhase=Plot[deltaArgand[m, mpi, mpi, Tmat[m*m + I 0.000001][[1,1]]],{m, mpi+mpi, 1.9}, AxesLabel->{M  [GeV/(c c)], \[Delta](\[Pi]\[Pi] -> \[Pi]\[Pi])}];
Export["f0PiPitoPiPiPhase.pdf", f0PiPitoPiPiPhase];

f0PiPitoPiPiArgand=ParametricPlot[{rho[m*m,mpi,mpi] Re[Tmat[m*m + I 0.000001][[1,1]]], rho[m*m,mpi,mpi] Im[Tmat[m*m + I 0.000001][[1,1]]]},{m, mpi+mpi, 1.9}, AxesLabel->{Re[T(\[Pi]\[Pi] -> \[Pi]\[Pi])], Im[T(\[Pi]\[Pi] -> \[Pi]\[Pi])]}];
Export["f0PiPitoPiArgand.pdf", f0PiPitoPiPiArgand];

f0PiPitoKKTsqr=Plot[ { rho[m*m,mpi,mpi] rho[m*m,mKp,mK0] Norm[Tmat[m*m + I 0.0000001][[1,3]]] Norm[Tmat[m*m + I 0.0000001][[1,3]]]}, {m, mKp+mK0, 1.9}, AxesLabel->{M [GeV/(c c)], "\[Rho](\[Pi]\[Pi])" "\[Rho](KK)" Abs[T(\[Pi]\[Pi] -> KK)] Abs[T(\[Pi]\[Pi] -> KK)]}];
Export["f0PiPitoKKArgandUnits.pdf", f0PiPitoKKTsqr];

f0PiPitoEtaEtaTsqr=Plot[ { rho[m*m,mpi,mpi] rho[m*m,meta,meta] Norm[Tmat[m*m + I 0.0000001][[1,4]]] Norm[Tmat[m*m + I 0.0000001][[1,4]]] }, {m, meta+meta, 1.9}, AxesLabel->{M  [GeV/(c c)], "\[Rho](\[Pi]\[Pi])" "\[Rho](\[Eta]\[Eta])" Abs[T(\[Pi]\[Pi] -> \[Eta]\[Eta])] Abs[T(\[Pi]\[Pi] -> \[Eta]\[Eta])]}];
Export["f0PiPitoEtaEtaArgandUnits.pdf", f0PiPitoEtaEtaTsqr];

f0PiPitoEtaEtaprimeTsqr=Plot[ { rho[m*m,mpi,mpi] rho[m*m,meta,metaprime] Norm[Tmat[m*m + I 0.0000001][[1,5]]] Norm[Tmat[m*m + I 0.0000001][[1,5]]] }, {m, meta+metaprime, 1.9}, AxesLabel->{M  [GeV/(c c)], "\[Rho](\[Pi]\[Pi])" "\[Rho](\[Eta]\[Eta]')" Abs[T(\[Pi]\[Pi] -> \[Eta]\[Eta]')] Abs[T(\[Pi]\[Pi] -> \[Eta]\[Eta]')]}];
Export["f0PiPitoEtaEtaprimeArgandUnits.pdf", f0PiPitoEtaEtaprimeTsqr];

