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
 This Mathematica-script calculates some relevant quantities for the scattering process based on the K-matrix parametrization of the pi1-wave published in Eur. Phys. J. C81 (2021) no.12, 1056; doi:10.1140/epjc/s10052-021-09821-2; arXiv:2008.11566.
 *)

Get[FileNameJoin[{Directory[], "KmatUtils.m"}]]

(* Precision[1.2]; *)
(* Precision[12/10]; *)

(* masses of the decay products *)
mpi:=0.13957;
meta:=0.547862;
metaprime:=0.95778;

(* K-matrix parameter *)
mpi11600:=1.385518356865105
gpieta:= 0.8056407725708118
gpietaprime:= 1.046952875320115
c00:= 1.0500
c01:= 0.1516321466620786
c11:= -0.2461137125335851

q0PiEta:=BreakupMomQDefaultFromS[mpi11600 mpi11600,mpi mpi, meta meta]; 
q0PiEtaprime:=BreakupMomQDefaultFromS[mpi11600 mpi11600,mpi mpi, metaprime metaprime]; 


KmatPi1BW[s_]:=
  Module[{resultMatr},
	 qPiEta = BreakupMomQDefaultFromS[s,mpi mpi, meta meta];
	 qPiEtaprime = BreakupMomQDefaultFromS[s,mpi mpi, metaprime metaprime];
	 BWPiEta = BWBarrier[1, qPiEta, q0PiEta];
	 BWPiEtaprime = BWBarrier[1, qPiEtaprime, q0PiEtaprime]; 
	 resultMatr = {{BWPiEta*(gpieta*gpieta/(mpi11600*mpi11600-s)+c00)*BWPiEta, BWPiEta*(gpieta*gpietaprime/(mpi11600*mpi11600-s)+c01)*BWPiEtaprime},
		       {BWPiEtaprime*(gpietaprime*gpieta/(mpi11600*mpi11600-s)+c01)*BWPiEta, BWPiEtaprime*(gpietaprime*gpietaprime/(mpi11600*mpi11600-s)+c11)*BWPiEtaprime}};
	 resultMatr
	 ];
         

ChewMmat[s_]:={{c[s,mpi,meta], 0.},{0., c[s,mpi,metaprime]}};

KChewMmat[s_]:=KmatPi1BW[s].ChewMmat[s];

Idmat:=IdentityMatrix[2];      
IKChewMmat[s_]:= Idmat + KChewMmat[s];

invIKChewMmat[s_]:= Inverse[IKChewMmat[s]];
	
Tmat[s_]:=invIKChewMmat[s].KmatPi1BW[s];

plotPiEtaPiEtaTreal=Plot[Re[Tmat[b*b + I 0.000001][[1,1]]],{b, mpi+meta, 2.0}, AxesLabel->{M [GeV/(c c)], Re[T(\[Pi]\[Eta] -> \[Pi]\[Eta])]}];
Export["pi1PiEtaPiEtaTreal.pdf", plotPiEtaPiEtaTreal];
plotPiEtaPiEtaTimag=Plot[Im[Tmat[b*b + I 0.000001][[1,1]]],{b, mpi+meta, 2.0}, AxesLabel->{M [GeV/(c c)], Im[T(\[Pi]\[Eta] -> \[Pi]\[Eta])]}];
Export["pi1PiEtaPiEtaTimag.pdf", plotPiEtaPiEtaTimag];

plotPiEtaPiEtaprimeTreal=Plot[Re[Tmat[b*b + I 0.000001][[2,1]]],{b, mpi+metaprime, 2.0}, AxesLabel->{M [GeV/(c c)], Re[T(\[Pi]\[Eta] -> \[Pi]\[Eta]')]}];
Export["pi1PiEtaPiEtaprimeTreal.pdf", plotPiEtaPiEtaprimeTreal];
plotPiEtaPiEtaprimeTimag=Plot[Im[Tmat[b*b + I 0.000001][[2,1]]],{b, mpi+metaprime, 2.0}, AxesLabel->{M [GeV/(c c)], Im[T(\[Pi]\[Eta] -> \[Pi]\[Eta]')]}];
Export["pi1PiEtaPiEtaprimeTimag.pdf", plotPiEtaPiEtaprimeTimag];

plotPiEtaprimePiEtaprimeTreal=Plot[Re[Tmat[b*b + I 0.000001][[2,2]]],{b, mpi+metaprime, 2.0}, AxesLabel->{M [GeV/(c c)], Re[T(\[Pi]\[Eta]' -> \[Pi]\[Eta]')]}];
Export["pi1PiEtaprimePiEtaprimeTreal.pdf", plotPiEtaprimePiEtaprimeTreal];
plotPiEtaprimePiEtaprimeTimag=Plot[Im[Tmat[b*b + I 0.000001][[2,2]]],{b, mpi+metaprime, 2.0}, AxesLabel->{M [GeV/(c c)], Im[T(\[Pi]\[Eta]' -> \[Pi]\[Eta]')]}];
Export["pi1PiEtaprimePiEtaprimeTimag.pdf", plotPiEtaprimePiEtaprimeTimag];


plotPiEtaPiEtaElasticity=Plot[Norm[Sii[b*b + I 0.0001, mpi, meta, Tmat[b*b + I 0.000001][[1,1]]]],{b, mpi+meta, 2.0}, AxesLabel->{M [GeV/(c c)], \[Eta](\[Pi]\[Eta] -> \[Pi]\[Eta])}];
Export["pi1PiEtaPiEtaElasticity.pdf", plotPiEtaPiEtaElasticity];

plotPiEtaPiEtaPhase=Plot[deltaArgand[m, mpi, meta, Tmat[m*m + I 0.000001][[1,1]]],{m, mpi+meta, 2.0}, AxesLabel->{M, \[Delta](\[Pi]\[Eta] -> \[Pi]\[Eta])}];
Export["pi1PiEtaPiEtaPhase.pdf", plotPiEtaPiEtaPhase];

plotPiEtaPiEtaArgand=ParametricPlot[{rho[m*m,mpi,meta] Re[Tmat[m*m + I 0.000001][[1,1]]], rho[m*m,mpi,meta] Im[Tmat[m*m + I 0.0001][[1,1]]]},{m, mpi+meta, 2.0}];
Export["pi1PiEtaPiEtaArgand.pdf", plotPiEtaPiEtaArgand];

plotPiEtaprimePiEtaprimeArgand=ParametricPlot[{rho[m*m,mpi,metaprime] Re[Tmat[m*m + I 0.000001][[2,2]]], rho[m*m,mpi,metaprime] Im[Tmat[m*m + I 0.0001][[2,2]]]},{m, mpi+metaprime, 2.0}];
Export["pi1PiEtaprimePiEtaprimeArgand.pdf", plotPiEtaprimePiEtaprimeArgand];

plotPiEtaprimePiEtaprimePhase=Plot[deltaArgand[m, mpi, metaprime, Tmat[m*m + I 0.000001][[2,2]]],{m, mpi+metaprime, 2.0}, AxesLabel->{M, \[Delta](\[Pi]\[Eta]' -> \[Pi]\[Eta]')}];
Export["pi1PiEtaprimePiEtaprimePhase.pdf", plotPiEtaprimePiEtaprimePhase];

plotPiEtaprimePiEtaprimeElasticity=Plot[Norm[Sii[b*b + I 0.000001, mpi, metaprime, Tmat[b*b + I 0.000001][[2,2]]]],{b, mpi+metaprime, 2.0}, AxesLabel->{M [GeV/(c c)], \[Eta](\[Pi]\[Eta]' -> \[Pi]\[Eta]')}];
Export["pi1PiEtaprimePiEtaprimeElasticity.pdf", plotPiEtaprimePiEtaprimeElasticity];
