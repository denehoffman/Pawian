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

***********************************************************************************************

    This Skript is meant to illustrate the shape of the stable and smeared Chew-Manedlstam  (CM) function in case of an unstable resonance. The Chew-Mandelstam function used here is based on Dudek et al  Phys. Rev. D91, 054008 (2015), Appendix B. The generalization towards unstable particles is based on Basdevant & Berger, Phys.Rev.D 19 (1979) 239. 
    The nomenclature is as follows:
    Imagine two particle scattering: P1* P2 -> P1* P2 
    where P1* is unstable and decays to p1, p2
    In the following Skript the example of rho-pi, rho-> pi+ pi- is implemented for illustration purposes..
    Otherwise, the user needs to choose the following input parameters:
    M1: Breit-Wigner/PDG mass of the unstable resonance
    m1,m2: Daughter particles of the decayin resonance
    fPole, mPole: Pole position of the unstable resonance in the complex energy plane as pbtained using the skript PoleExtract_CM_Dudek.m (see header for more information)
*)

sqrts := x + I y;
y=0;
s := sqrts * sqrts;

(* Rho-Pi Case is implemented for demonstration, if needed change input values here *)
mpi0:=0.1349768;
mpi:=0.13957;

(* Rho Mass of the unstable resonance*)
M1:=0.77526;

(* Masses of daughters *)
m1:=mpi;
m2:=mpi;

(* Pole Position as obtained from PoleExtract_CM_Dudek.m using the values for rho->pi+ pi-. The calculation in PoleExtract_CM_Dudek.m relies on finding the singularity in ds. ds is the denominator of dispersion relation, thus its singularity defines the defines pole position, see Basdevant & Berger, Phys.Rev.D 19 (1979) 239  eq(12) and eq(14)*)
fPole:=Sqrt[0.3529];
mPole:=Sqrt[0.8161];

threshold:=(m1+m2)^2;

(* Definition of the CM function after Dudek et al  Phys. Rev. D91, 054008 (2015), sprime and scaler act as coordinate transformations to shift the integral range from [threshold, infinity] to [0, 1] which saves computation time *)
rho[s_,m1_,m2_]:=Sqrt[1-(m1+m2)^2/s] Sqrt[1-(m1-m2)^2/s];
xi[s_,m1_,m2_]:=1-(m1+m2)^2/s;
sprime[a_,t_]:=(a+(t/(1-t)));
scaler[u_]:=(1/((1-u)^2));

c[s_,m_,u_]:=rho[s,m,u]/\[Pi] Log[(xi[s,m,u]+rho[s,m,u])/(xi[s,m,u]-rho[s,m,u])]-xi[s,m,u]/\[Pi] (u-m)/(u+m) Log[u/m];

(* taking care of singularities and numerical stability *)
cm[S_,m_,u_]:= Piecewise[{{c[S,m,u],S!=0&&1-(m+u)^2/S-Sqrt[1-(m-u)^2/S] Sqrt[1-(m+u)^2/S]!=0&&(1-(m+u)^2/S+Sqrt[1-(m-u)^2/S] Sqrt[1-(m+u)^2/S])/(1-(m+u)^2/S-Sqrt[1-(m-u)^2/S] Sqrt[1-(m+u)^2/S])!=0},{0,True}}];

ds[t_,m1_,m2_,fR_,mR_,a_]:=sprime[a,t]-mR^2+(fR^2*(sprime[a,t]-a)*cm[sprime[a,t],m1,m2]);
Sigma[a_,t_,m1_,m2_]:=(sprime[a,t]-a)*cm[sprime[a,t],m1,m2];

(* change PrecisionGoal if needed, e-4 is a sufficient compromise for plot purposes (although this can sstill take several minutes) *)
Ctilde[S_,mPole_,fPole_,u_,m1_,m2_]:=NIntegrate[SetPrecision[-(((1/\[Pi])*fPole^2*Im[Sigma[threshold,s1,m1,m2]])/Abs[ds[s1,m1,m2,fPole,mPole,threshold]]^2)*cm[S,Sqrt[sprime[threshold,s1]],u]*scaler[s1],30],{s1,0,1},Method->{"GlobalAdaptive","SingularityHandler"->"DoubleExponential"},MaxRecursion->20,MaxPoints->10000000,WorkingPrecision->25,PrecisionGoal->4];

(* Finally plot and save it*)
pl2=ReImPlot[c[s,M1,mpi0],{x,0.01,3},PlotLegends->Automatic, ReImStyle->{ColorData[97,3],ColorData[97,7]},AxesLabel->{Sqrt[S] [GeV],C}];


pl=ReImPlot[-Ctilde[s,mPole,fPole,mpi0,m1,m2],{x,0.01,3},MaxRecursion->0,PlotLegends->Automatic, ReImStyle->{ColorData[97,2],ColorData[97,14]},AxesLabel->{Sqrt[S][GeV],OverTilde[C]}];

(* Comment the following lines if you do not want to save the plots *)
Export["Stable_CM_Dudek.pdf",pl2];
Export["Unstable_CM_Dudek.pdf",pl];
pl3=Show[pl2,pl];
Export["CM_Dudek_compare.pdf",pl3];
