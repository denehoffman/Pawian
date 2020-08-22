(*  Copyright 2020 
    Meike Kuessner
    mail: mkuessner@ep1.rub.de
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

***********************************************************************
This Script computes the Pole Position of a resonance by its mass,
width and masses of the decay products in case of a two-body decay. 
The calculation is based on Phys.Rev.D 19 (1979) 239  eq(12) and eq(14) 
and the Chew-Mandelstam function is based on 
Dudek et al  Phys. Rev. D91, 054008 (2015), Appendix B
***********************************************************************
*)

sqrts := x + I y;
y=0;
s := sqrts sqrts ;

(* Rho-Pi Case is implemented for demonstration, if needed change input values here *)

(* Mothers Mass and Width *)
mpi:=0.13957;
M1:=0.77526;
\[CapitalGamma]:=0.1491; (*Width of unstable resonance*)
(* Masses of daughters *)
m1:=mpi; 
m2:=mpi;

(* Define CM Function based on Dudek et al  Phys. Rev. D91, 054008 (2015), Appendix B *)
rho[s_,m1_,m2_]:=Sqrt[1-(m1+m2)^2/s] Sqrt[1-(m1-m2)^2/s];
xi[s_,m1_,m2_]:=1-(m1+m2)^2/s;

c[s_,m_,u_]:=rho[s,m,u]/\[Pi] Log[(xi[s,m,u]+rho[s,m,u])/(xi[s,m,u]-rho[s,m,u])]-xi[s,m,u]/\[Pi] (u-m)/(u+m) Log[u/m];

(* Singularity of denominator |d(s)|^2 of dispersion relation defines pole position, see Basdevant & Berger, Phys.Rev.D 19 (1979) 239  eq(12) and eq(14) *)
ds2[s_,m1_,m2_,fR_,mR_]:=Abs[s-mR+fR*(s-(m1+m2)^2)*c[s,m1,m2]]^2;

(* Determine Minimum and save *)
result=FindMinimum[SetPrecision[ds2[(M1-I*\[CapitalGamma]/2)^2,m2,m2,a,b],30],{a,b},WorkingPrecision->20];
{fsq,msq}={a,b}/.result[[2]];
Print["m^2 = ", msq, " ", "f^2 = ", fsq];
Print["m = ", Sqrt[msq], " ", "f = ", Sqrt[fsq]];

(* Plot Contour around Minimum and Save Plot, change plot style here if needed or comment to not save the plot*)
plot=ContourPlot[Log[ds2[(M1-I*\[CapitalGamma]/2)^2,m2,m2,f,m]],{f,fsq-0.1,fsq+0.1},{m,msq-0.1,msq+0.1},Contours->100,PlotLegends->Automatic,WorkingPrecision->30,ContourStyle->None, ColorFunction->"Rainbow", Epilog->({Black,PointSize[.03],Point[{fsq,msq}],Text["minimum",{fsq,msq*0.975}]}),FrameLabel->{f^2,m^2}];

Export["Pole.pdf",plot];
