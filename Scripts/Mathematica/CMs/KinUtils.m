(*  Copyright 2025                                                                                
 Bertram Kopf (bertram@ep1.rub.de)
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
  This Mathematica-script is a colletion of functions which are needed to perform calculations with the K-matrix formalism.
 *)

(* phase space factor *)
rho[s_,m1_,m2_]:= Sqrt[1-(m1+m2)^2/s] Sqrt[1-(m1-m2)^2/s];

(* break up momentum *)
BreakupMomQDefaultFromS[s_, s1_, s2_]:=
  Module[{result},
	 result=rho[s, Sqrt[s1], Sqrt[s2]] Sqrt[s]/2.;
	 (* result = Sqrt[s*s + s1*s1 + s2*s2 - 2. * ( s*s1 + s*s2 + s1*s2)]/Sqrt[4.*s]; *)
	 result
	 ];

(* CM Function based on Dudek et al  Phys. Rev. D91, 054008 (2015), Appendix B; subtraction constant at s=s_thresh CM(s_thresh)=0 *)
xi[s_,m1_,m2_]:=1-(m1+m2)^2/s;
cmDudekPre[s_,m_,u_] := rho[s,m,u]/\[Pi] Log[(xi[s,m,u]+rho[s,m,u])/(xi[s,m,u]-rho[s,m,u])]-xi[s,m,u]/\[Pi] (u-m)/(u+m) Log[u/m];

(* making sure to be in the 1 quadrant for s=real *)
cmDudek[s_,m_,u_] := cmDudekPre[s + I * If[Im[s] == 0, 10.^-8, 0.], m, u];


(* CM Function based on Reid with subtraction constant at s=0 CM(0)=0 *)
(* cmReid[s_,m_,u_]:= - ( cmDudek[s,m,u] - 1/\[Pi] (1. + 2 m u / (u u - m m) Log[u/m] ) ); *)

(* proteted case for u==m with L'Hospital *) 
cmReid[s_,m_,u_]:= - ( cmDudek[s,m,u] - If[m != u, 1/\[Pi] (1. + 2 m u / (u u - m m) Log[u/m]), 2./\[Pi]] );

(* cmReid[s_,m_,u_]:= - ( cmDudek[s,m,u] - If[m != u, 1/\[Pi] (1. + 2 m u / (u u - m m) Log[u/m] ), 2./\[Pi]  *)


(* Blatt-Weisskopf factors *)
BlattWeisskopf[l_, z_]:=
  Module[{result},
	 result = Which[l==0, 1.0 + 0. I, l==1, Sqrt[2. z / (z+1.)] + 0. I, l==2, Sqrt[13. z z / (z z + 3. z + 9.)]+ 0. I];
	 result
	 ];

BWBarrier[l_, q_, q0_]:=
  Module[{result},
	 qR:=0.1973;
  z:=(Norm[q] Norm[q])/(Norm[qR] Norm[qR]);
	 (* Print["z: ", z]; *) 
	 z0:=(Norm[q0] Norm[q0])/(Norm[qR] Norm[qR]);
	 (* Print["z0: ", z0]; *)
	 result:= BlattWeisskopf[l,z]/BlattWeisskopf[l,z0];
	 result
	 ];

(* elasticity *)
Sii[s_, m1_, m2_, Tii_] := 1. + 2. I Sqrt[Re[rho[s,m1,m2]]] Tii Sqrt[Re[rho[s,m1,m2]]];


(* extract delta phase of an Argand plot *)
deltaArgand[mass_, m1_, m2_, Tii_]:=
  Module[{delta},
	 Tiirho :=rho[mass*mass,m1,m2] Tii;
  ReE := Re[Tiirho];
  ImE := Im[Tiirho] - 0.5;
	 tmpdelta := ArcTan[Abs[ReE], ImE];
	 delta := 0.5 tmpdelta 180./Pi + 45.0;
	 Which[ReE  < 0.0, delta = 180.0 - delta];
	 delta
	 ];

(* Breit-Wigner function w/ Blatt-Weisskopf factors *)
BreitWignerBlattWRel[sCurrent_, orbMom_, M0_, Width_, Ma_, Mb_]:=
  Module[{BW},
	 s0=M0 M0;
	 rho0=rho[s0, Ma, Mb];
	 rhoCurrent=rho[sCurrent, Ma, Mb];
         momQ0=BreakupMomQDefaultFromS[s0, Ma Ma, Mb Mb];
	 momQCurrent=BreakupMomQDefaultFromS[sCurrent, Ma Ma, Mb Mb];
	 bwBarrierFac=If[orbMom==0, 1., BWBarrier[orbMom, momQCurrent, momQ0]];
	 BW = M0 Width  bwBarrierFac / (s0 - sCurrent - I (rhoCurrent/rho0) M0 Width bwBarrierFac bwBarrierFac); 
	 BW
	 ];

BreitWignerBlattWRelNorm[sCurrent_, orbMom_, M0_, Width_, Ma_, Mb_]:=
  Module[{BW},
         s0=M0 M0;
         rho0=rho[s0, Ma, Mb];
         rhoCurrent=rho[sCurrent, Ma, Mb];
         momQ0=BreakupMomQDefaultFromS[s0, Ma Ma, Mb Mb];
         momQCurrent=BreakupMomQDefaultFromS[sCurrent, Ma Ma, Mb Mb];
         bwBarrierFac=If[orbMom==0, 1., BWBarrier[orbMom, momQCurrent, momQ0]];
         GammaM= Width (rhoCurrent/rho0) bwBarrierFac bwBarrierFac;
	 BW = Sqrt[M0 GammaM] / (s0 - sCurrent - I M0 GammaM);
         BW
         ];

(* sinplified barrier function based on s *)
BarrierSimple[s_, sthresh_, s0_, orbMom_]:= If[orbMom==0, 1., ( (s-sthresh) / (s+s0) )^orbMom ];

