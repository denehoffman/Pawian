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
  This Mathematica-script is a colletion of functions which are needed to perform calculations with the K-matrix formalism.
 *)

(* phase space factor *)
rho[s_,m1_,m2_]:=Sqrt[1-(m1+m2)^2/s] Sqrt[1-(m1-m2)^2/s];

(* break up momentum *)
BreakupMomQDefaultFromS[s_, s1_, s2_]:=
  Module[{result},
	 result=rho[s, Sqrt[s1], Sqrt[s2]] Sqrt[s]/2.;
	 (* result = Sqrt[s*s + s1*s1 + s2*s2 - 2. * ( s*s1 + s*s2 + s1*s2)]/Sqrt[4.*s]; *)
	 result
	 ];

(* CM Function based on Dudek et al  Phys. Rev. D91, 054008 (2015), Appendix B *)
xi[s_,m1_,m2_]:=1-(m1+m2)^2/s;
c[s_,m_,u_]:=rho[s,m,u]/\[Pi] Log[(xi[s,m,u]+rho[s,m,u])/(xi[s,m,u]-rho[s,m,u])]-xi[s,m,u]/\[Pi] (u-m)/(u+m) Log[u/m];

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
Sii[s_, m1_, m2_, Tii_] := 1. + 2 I Sqrt[Re[rho[s,m1,m2]]] Tii Sqrt[Re[rho[s,m1,m2]]];


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
