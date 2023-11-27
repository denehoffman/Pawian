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
 This Mathematica-script calculates some relevant quantities for the scattering process based on the K-matrix parametrization of the f2-wave published in Eur. Phys. J. C81 (2021) no.12, 1056; doi:10.1140/epjc/s10052-021-09821-2; arXiv:2008.11566.
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

(* K-matrix parameter *)
mf21270 := 1.152989902015372;
mf21525 := 1.483591387238627;
mf21810 := 1.729226406267156;
mf22100 := 1.969999999998483; 

gpipif21270 := 0.4003279488527079;
g4pif21270  := 0.1547855579730676;
gKKf21270 := -0.08900000050965565;
getaetaf21270 := -0.001134042567039056;
gpipif21525 := 0.01819999860301777;
g4pif21525 := 0.1729999857292003;
gKKf21525 := 0.3239256460800748;
getaetaf21525 := 0.1525632940358444;
gpipif21810 := -0.06709734248406923;
g4pif21810 := 0.2294058495107097;
gKKf21810 := -0.4313312049403804;
getaetaf21810 := 0.2372109717029494;
gpipif22100 := -0.4992419169695315;
g4pif22100 := 0.1929461301610656;
gKKf22100 := 0.279745459722337;
getaetaf22100 := -0.03986921341499428;

c00 := -0.04318879859335506;
c01 :=   0.;
c02 :=   0.00984386157861733;
c03 :=   0.01027674144082516;
c11 :=   0.;
c12 :=   0.;
c13 :=   0.;
c22 :=   -0.07343513098907345;
c23 :=   0.05532677058094471;
c33 :=   -0.05183689223881875;

q0PiPif21270:=BreakupMomQDefaultFromS[mf21270 mf21270,mpi mpi, mpi mpi];
q04Pif21270:=BreakupMomQDefaultFromS[mf21270 mf21270,m2pi m2pi, m2pi m2pi];
q0KKf21270:=BreakupMomQDefaultFromS[mf21270 mf21270,mKp mKp, mK0 mK0];
q0EtaEtaf21270:=BreakupMomQDefaultFromS[mf21270 mf21270,meta meta, meta meta];

q0PiPif21525:=BreakupMomQDefaultFromS[mf21525 mf21525,mpi mpi, mpi mpi];
q04Pif21525:=BreakupMomQDefaultFromS[mf21525 mf21525,m2pi m2pi, m2pi m2pi];
q0KKf21525:=BreakupMomQDefaultFromS[mf21525 mf21525,mKp mKp, mK0 mK0];
q0EtaEtaf21525:=BreakupMomQDefaultFromS[mf21525 mf21525,meta meta, meta meta];

q0PiPif21810:=BreakupMomQDefaultFromS[mf21810 mf21810,mpi mpi, mpi mpi];
q04Pif21810:=BreakupMomQDefaultFromS[mf21810 mf21810,m2pi m2pi, m2pi m2pi];
q0KKf21810:=BreakupMomQDefaultFromS[mf21810 mf21810,mKp mKp, mK0 mK0];
q0EtaEtaf21810:=BreakupMomQDefaultFromS[mf21810 mf21810,meta meta, meta meta];

q0PiPif22100:=BreakupMomQDefaultFromS[mf22100 mf22100,mpi mpi, mpi mpi];
q04Pif22100:=BreakupMomQDefaultFromS[mf22100 mf22100,m2pi m2pi, m2pi m2pi];
q0KKf22100:=BreakupMomQDefaultFromS[mf22100 mf22100,mKp mKp, mK0 mK0];
q0EtaEtaf22100:=BreakupMomQDefaultFromS[mf22100 mf22100,meta meta, meta meta];

Kmatf2BW[s_]:=
  Module[{resultMatr},
	 qPiPi = BreakupMomQDefaultFromS[Re[s], mpi mpi, mpi mpi];
	 q4Pi = BreakupMomQDefaultFromS[Re[s], m2pi m2pi, m2pi m2pi];
	 qKK = BreakupMomQDefaultFromS[Re[s], mKp mKp, mK0 mK0];
	 qEtaEta = BreakupMomQDefaultFromS[Re[s], meta meta, meta meta];
	 	 
	 BWPiPif21270 = BWBarrier[2, qPiPi, q0PiPif21270];
	 BW4Pif21270 = BWBarrier[2, q4Pi, q04Pif21270];
	 BWKKf21270 = BWBarrier[2, qKK, q0KKf21270];
	 BWEtaEtaf21270 = BWBarrier[2, qEtaEta, q0EtaEtaf21270];

	 BWPiPif21525 = BWBarrier[2, qPiPi, q0PiPif21525];
	 BW4Pif21525 = BWBarrier[2, q4Pi, q04Pif21525];
	 BWKKf21525 = BWBarrier[2, qKK, q0KKf21525];
	 BWEtaEtaf21525 = BWBarrier[2, qEtaEta, q0EtaEtaf21525];

	 BWPiPif21810 = BWBarrier[2, qPiPi, q0PiPif21810];
	 BW4Pif21810 = BWBarrier[2, q4Pi, q04Pif21810];
	 BWKKf21810 = BWBarrier[2, qKK, q0KKf21810];
	 BWEtaEtaf21810 = BWBarrier[2, qEtaEta, q0EtaEtaf21810];

	 BWPiPif22100 = BWBarrier[2, qPiPi, q0PiPif22100];
	 BW4Pif22100 = BWBarrier[2, q4Pi, q04Pif22100];
	 BWKKf22100 = BWBarrier[2, qKK, q0KKf22100];
	 BWEtaEtaf22100 = BWBarrier[2, qEtaEta, q0EtaEtaf22100];
	 
	 KPiPitoPiPi=BWPiPif21270 (gpipif21270 gpipif21270/(mf21270 mf21270-s)+c00) BWPiPif21270
	 +BWPiPif21525 (gpipif21525 gpipif21525/(mf21525 mf21525-s)+c00) BWPiPif21525
	 +BWPiPif21810 (gpipif21810 gpipif21810/(mf21810 mf21810-s)+c00) BWPiPif21810
	 +BWPiPif22100 (gpipif22100 gpipif22100/(mf22100 mf22100-s)+c00) BWPiPif22100;

	 KPiPito4Pi=BWPiPif21270 (gpipif21270 g4pif21270/(mf21270 mf21270-s)+c01) BW4Pif21270
	 +BWPiPif21525 (gpipif21525 g4pif21525/(mf21525 mf21525-s)+c01) BW4Pif21525
	 +BWPiPif21810 (gpipif21810 g4pif21810/(mf21810 mf21810-s)+c01) BW4Pif21810
	 +BWPiPif22100 (gpipif22100 g4pif22100/(mf22100 mf22100-s)+c01) BW4Pif22100;

	 KPiPitoKK=BWPiPif21270 (gpipif21270 gKKf21270/(mf21270 mf21270-s)+c02) BWKKf21270
	 +BWPiPif21525 (gpipif21525 gKKf21525/(mf21525 mf21525-s)+c02) BWKKf21525
	 +BWPiPif21810 (gpipif21810 gKKf21810/(mf21810 mf21810-s)+c02) BWKKf21810
	 +BWPiPif22100 (gpipif22100 gKKf22100/(mf22100 mf22100-s)+c02) BWKKf22100;

	 KPiPitoEtaEta=BWPiPif21270 (gpipif21270 getaetaf21270/(mf21270 mf21270-s)+c03) BWEtaEtaf21270
	 +BWPiPif21525 (gpipif21525 getaetaf21525/(mf21525 mf21525-s)+c03) BWEtaEtaf21525
	 +BWPiPif21810 (gpipif21810 getaetaf21810/(mf21810 mf21810-s)+c03) BWEtaEtaf21810
	 +BWPiPif22100 (gpipif22100 getaetaf22100/(mf22100 mf22100-s)+c03) BWEtaEtaf22100;

	 K4Pito4Pi=BW4Pif21270 (g4pif21270 g4pif21270/(mf21270 mf21270-s)+c11) BW4Pif21270
	 +BW4Pif21525 (g4pif21525 g4pif21525/(mf21525 mf21525-s)+c11) BW4Pif21525
	 +BW4Pif21810 (g4pif21810 g4pif21810/(mf21810 mf21810-s)+c11) BW4Pif21810
	 +BW4Pif22100 (g4pif22100 g4pif22100/(mf22100 mf22100-s)+c11) BW4Pif22100;

	 K4PitoKK=BW4Pif21270 (g4pif21270 gKKf21270/(mf21270 mf21270-s)+c12) BWKKf21270
	 +BW4Pif21525 (g4pif21525 gKKf21525/(mf21525 mf21525-s)+c12) BWKKf21525
	 +BW4Pif21810 (g4pif21810 gKKf21810/(mf21810 mf21810-s)+c12) BWKKf21810
	 +BW4Pif22100 (g4pif22100 gKKf22100/(mf22100 mf22100-s)+c12) BWKKf22100;

	 K4PitoEtaEta=BW4Pif21270 (g4pif21270 getaetaf21270/(mf21270 mf21270-s)+c13) BWEtaEtaf21270
	 +BW4Pif21525 (g4pif21525 getaetaf21525/(mf21525 mf21525-s)+c13) BWEtaEtaf21525
	 +BW4Pif21810 (g4pif21810 getaetaf21810/(mf21810 mf21810-s)+c13) BWEtaEtaf21810
	 +BW4Pif22100 (g4pif22100 getaetaf22100/(mf22100 mf22100-s)+c13) BWEtaEtaf22100;

	 KKKtoKK=BWKKf21270 (gKKf21270 gKKf21270/(mf21270 mf21270-s)+c22) BWKKf21270
	 +BWKKf21525 (gKKf21525 gKKf21525/(mf21525 mf21525-s)+c22) BWKKf21525
	 +BWKKf21810 (gKKf21810 gKKf21810/(mf21810 mf21810-s)+c22) BWKKf21810
	 +BWKKf22100 (gKKf22100 gKKf22100/(mf22100 mf22100-s)+c22) BWKKf22100;

	 KKKtoEtaEta=BWKKf21270 (gKKf21270 getaetaf21270/(mf21270 mf21270-s)+c23) BWEtaEtaf21270
	 +BWKKf21525 (gKKf21525 getaetaf21525/(mf21525 mf21525-s)+c23) BWEtaEtaf21525
	 +BWKKf21810 (gKKf21810 getaetaf21810/(mf21810 mf21810-s)+c23) BWEtaEtaf21810
	 +BWKKf22100 (gKKf22100 getaetaf22100/(mf22100 mf22100-s)+c23) BWEtaEtaf22100;

	 KEtaEtatoEtaEta=BWEtaEtaf21270 (getaetaf21270 getaetaf21270/(mf21270 mf21270-s)+c33) BWEtaEtaf21270
	 +BWEtaEtaf21525 (getaetaf21525 getaetaf21525/(mf21525 mf21525-s)+c33) BWEtaEtaf21525
	 +BWEtaEtaf21810 (getaetaf21810 getaetaf21810/(mf21810 mf21810-s)+c33) BWEtaEtaf21810
	 +BWEtaEtaf22100 (getaetaf22100 getaetaf22100/(mf22100 mf22100-s)+c33) BWEtaEtaf22100;
	 
	 resultMatr = {
	   {KPiPitoPiPi,   KPiPito4Pi,   KPiPitoKK,   KPiPitoEtaEta}, 
	   {KPiPito4Pi,    K4Pito4Pi,    K4PitoKK,    K4PitoEtaEta},
	   {KPiPitoKK,     K4PitoKK,     KKKtoKK,     KKKtoEtaEta},
	   {KPiPitoEtaEta, K4PitoEtaEta, KKKtoEtaEta, KEtaEtatoEtaEta}}; 
	   resultMatr
	   ]; 	  
         

	 ChewMmat[s_]:={{c[s,mpi,mpi], 0., 0., 0.},
			{0., c[s,m2pi,m2pi], 0. ,0.},
			{0., 0., c[s,mKp,mK0], 0.},
			{0., 0., 0., c[s,meta,meta]}
	 };


KChewMmat[s_]:=Kmatf2BW[s].ChewMmat[s];
(* Print["CM[1.8]: ", ChewMmat[1.8] // MatrixForm]; *)
      
      
Idmat:=IdentityMatrix[4];      
IKChewMmat[s_]:= Idmat + KChewMmat[s];

invIKChewMmat[s_]:= Inverse[IKChewMmat[s]];
	
Tmat[s_]:=invIKChewMmat[s].Kmatf2BW[s];

	 (*	 Print["Tmat[1.8]: ", Tmat[1.8] // MatrixForm]; *)  
		       
f2PiPitoPiPiTreal=Plot[Re[Tmat[b*b + I 0.000001][[1,1]]],{b, mpi+mpi, 1.9}, AxesLabel->{M [GeV/(c c)], Re[T(\[Pi]\[Pi] -> \[Pi]\[Pi])]}];
Export["f2PiPitoPiPiTreal.pdf", f2PiPitoPiPiTreal];
f2PiPitoPiPiTimag=Plot[Im[Tmat[b*b + I 0.000001][[1,1]]],{b, mpi+mpi, 1.9}, AxesLabel->{M [GeV/(c c)], Im[T(\[Pi]\[Pi] -> \[Pi]\[Pi])]}];
Export["f2PiPitoPiPiTimag.pdf", f2PiPitoPiPiTimag];

f2PiPitoEtaEtaTreal=Plot[Re[Tmat[b*b + I 0.000001][[1,4]]],{b, meta+meta, 1.9}, AxesLabel->{M [GeV/(c c)], ReT[(\[Pi]\[Pi] -> \[Eta]\[Eta])]}];
Export["f2PiPitoEtaEtaTreal.pdf", f2PiPitoEtaEtaTreal];
f2PiPitoEtaEtaTimag=Plot[Im[Tmat[b*b + I 0.000001][[1,4]]],{b, meta+meta, 1.9}, AxesLabel->{M [GeV/(c c)], Im[T(\[Pi]\[Pi] -> \[Eta]\[Eta])]}, WorkingPrecision -> 32, AccuracyGoal -> 16];
Export["f2PiPitoEtaEtaTimag.pdf", f2PiPitoEtaEtaTimag];

f2EtaEtatoEtaEtaTreal=Plot[Re[Tmat[b*b + I 0.000001][[4,4]]],{b, meta+meta, 1.9}, AxesLabel->{M [GeV/(c c)], Re[T(\[Eta]\[Eta] -> \[Eta]\[Eta])]}];
Export["f2EtaEtatoEtaEtaTreal.pdf", f2EtaEtatoEtaEtaTreal];
f2EtaEtatoEtaEtaTimag=Plot[Im[Tmat[b*b + I 0.000001][[4,4]]],{b, meta+meta, 1.9}, AxesLabel->{M [GeV/(c c)], Im[T(\[Eta]\[Eta] -> \[Eta]\[Eta])]}, WorkingPrecision -> 32, AccuracyGoal -> 16];
Export["f2EtaEtatoEtaEtaTimag.pdf", f2EtaEtatoEtaEtaTimag];


f2PiPiElasticity=Plot[Norm[Sii[b*b + I 0.000001, mpi, mpi, Tmat[b*b + I 0.000001][[1,1]]]],{b, mpi+mpi, 1.9}, AxesLabel->{M [GeV/(c c)], \[Eta](\[Pi]\[Pi] -> \[Pi]\[Pi])}];
Export["f2PiPiElasticity.pdf", f2PiPiElasticity];

f2PiPitoPiPiPhase=Plot[deltaArgand[m, mpi, mpi, Tmat[m*m + I 0.000001][[1,1]]],{m, mpi+mpi, 1.9}, AxesLabel->{M  [GeV/(c c)], \[Delta](\[Pi]\[Pi] -> \[Pi]\[Pi])}];
Export["f2PiPitoPiPiPhase.pdf", f2PiPitoPiPiPhase];

f2PiPitoPiPiArgand=ParametricPlot[{rho[m*m,mpi,mpi] Re[Tmat[m*m + I 0.000001][[1,1]]], rho[m*m,mpi,mpi] Im[Tmat[m*m + I 0.000001][[1,1]]]},{m, mpi+mpi, 1.9}, AxesLabel->{Re[T(\[Pi]\[Pi] -> \[Pi]\[Pi])], Im[T(\[Pi]\[Pi] -> \[Pi]\[Pi])]}];
Export["f2PiPitoPiArgand.pdf", f2PiPitoPiPiArgand];

f2PiPitoKKTsqr=Plot[ { 5. rho[m*m,mpi,mpi] rho[m*m,mKp,mK0] Norm[Tmat[m*m + I 0.0000001][[1,3]]] Norm[Tmat[m*m + I 0.0000001][[1,3]]]}, {m, mKp+mK0, 1.9}, AxesLabel->{M [GeV/(c c)], 5 "\[Rho](\[Pi]\[Pi])" "\[Rho](KK)" Abs[T(\[Pi]\[Pi] -> KK)] Abs[T(\[Pi]\[Pi] -> KK)]}];
Export["f2PiPitoKKArgandUnits.pdf", f2PiPitoKKTsqr];

f2PiPitoEtaEtaTsqr=Plot[ { 5. rho[m*m,mpi,mpi] rho[m*m,meta,meta] Norm[Tmat[m*m + I 0.0000001][[1,4]]] Norm[Tmat[m*m + I 0.0000001][[1,4]]] }, {m, meta+meta, 1.9}, AxesLabel->{M  [GeV/(c c)], 5 "\[Rho](\[Pi]\[Pi])" "\[Rho](\[Eta]\[Eta])" Abs[T(\[Pi]\[Pi] -> \[Eta]\[Eta])] Abs[T(\[Pi]\[Pi] -> \[Eta]\[Eta])]}];
Export["f2PiPitoEtaEtaArgandUnits.pdf", f2PiPitoEtaEtaTsqr];

