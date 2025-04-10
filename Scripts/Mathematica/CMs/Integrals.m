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

Get[FileNameJoin[{Directory[], "../KinUtils.m"}]]

(* Precision[1.2]; *)
(* Precision[12/10]; *)

dsDudek[sprime_,m1_,m2_,f2R_,m2R_, s0_, orbL_] := m2R - sprime + f2R BarrierSimple[sprime, (m1+m2)^2, s0, orbL] cmDudek[sprime,m1,m2];

ImagSigmads2Dudek[sprime_,m1_,m2_,f2R_,m2R_, s0_, orbL_] := If[Abs[dsDudek[sprime,m1,m2,f2R,m2R, s0, orbL]]^2 < 10^30, - BarrierSimple[sprime, (m1+m2)^2, s0, orbL] Im[cmDudek[sprime, m1, m2]] / Abs[dsDudek[sprime,m1,m2,f2R,m2R, s0, orbL]]^2, 0.];

RNIntDudek[m1_, m2_, f2R_, m2R_, s0_, orbL_, scalingFac_, smax_] :=
  scalingFac /\[Pi] NIntegrate[ f2R ImagSigmads2Dudek[sb,m1,m2,f2R,m2R, s0, orbL]
			, {sb, (m1 + m2)^2, smax}, 
			Method -> {"GaussKronrodRule", "SymbolicProcessing" -> 0},
			MaxRecursion -> 20, 
			AccuracyGoal -> 5];


dsBB[sprime_,m1_,m2_,f2R_,m2R_, s0_, orbL_] := m2R - sprime - f2R BarrierSimple[sprime, (m1+m2)^2, s0, orbL] cmReid[sprime,m1,m2];

ImagSigmads2BB[sprime_,m1_,m2_,f2R_,m2R_, s0_, orbL_] := If[Abs[dsBB[sprime,m1,m2,f2R,m2R, s0, orbL]]^2 < 10^30, BarrierSimple[sprime, (m1+m2)^2, s0, orbL] Im[cmReid[sprime, m1, m2]] / Abs[dsBB[sprime,m1,m2,f2R,m2R, s0, orbL]]^2, 0.];

RNIntBB[m1_, m2_, f2R_, m2R_, s0_, orbL_, scalingFac_, smax_] :=
  scalingFac /\[Pi] NIntegrate[ f2R ImagSigmads2BB[sb,m1,m2,f2R, m2R, s0, orbL]
			, {sb, (m1 + m2)^2, smax}, 
			Method -> {"GaussKronrodRule", "SymbolicProcessing" -> 0},
			MaxRecursion -> 20, 
			AccuracyGoal -> 5];

RNIntBW[M0_, Width_, m1_, m2_, orbMom_, scalingFac_, smax_] :=
  scalingFac /\[Pi] NIntegrate[  rho[sb, m1, m2] Abs[BreitWignerBlattWRelNorm[sb, orbMom, M0, Width, m1, m2]]^2
                  , {sb, (m1 + m2)^2, smax}, 
                        Method -> {"GaussKronrodRule", "SymbolicProcessing" -> 0},
                        MaxRecursion -> 20, 
                        AccuracyGoal -> 5]; 

RNIntDudekviaBB[s_][m1_, m2_, mRecoil_, f2R_, m2R_, s0_, orbL_, scalingFac_] :=
  scalingFac /\[Pi] NIntegrate[ f2R ImagSigmads2BB[sb,m1,m2,f2R, m2R, s0, orbL] cmDudek[s,Sqrt[sb], mRecoil]
                        , {sb, (m1 + m2)^2, Infinity}, 
                        Method -> {"GaussKronrodRule", "SymbolicProcessing" -> 0},
                        MaxRecursion -> 20, 
                        AccuracyGoal -> 5];
 
RNIntDudekviaDudek[s_][m1_, m2_, mRecoil_, f2R_, m2R_, s0_, orbL_, scalingFac_] :=
  scalingFac /\[Pi] NIntegrate[ f2R ImagSigmads2Dudek[sb,m1,m2,f2R, m2R, s0, orbL] cmDudek[s,Sqrt[sb], mRecoil]
                        , {sb, (m1 + m2)^2, Infinity}, 
                        Method -> {"GaussKronrodRule", "SymbolicProcessing" -> 0},
                        MaxRecursion -> 20, 
                        AccuracyGoal -> 5];

RNIntDudekviaBW[s_][m1_, m2_, mrec_, orbMom_, M0_, Width_, scale_] :=
  scale /\[Pi] NIntegrate[ rho[sb, m1, m2] Abs[BreitWignerBlattWRelNorm[sb, orbMom, M0, Width, m1, m2]]^2 cmDudek[s, Sqrt[sb], mrec]
                        , {sb, (m1 + m2)^2, Infinity}, 
                        Method -> {"GaussKronrodRule", "SymbolicProcessing" -> 0},
                        MaxRecursion -> 20, 
                        AccuracyGoal -> 5];

RNIntBBviaBB[s_][m1_, m2_, mRecoil_, f2R_, m2R_, s0_, orbL_, scalingFac_] :=
  scalingFac /\[Pi] NIntegrate[ f2R ImagSigmads2BB[sb,m1,m2,f2R, m2R, s0, orbL] cmReid[s,Sqrt[sb], mRecoil]
                        , {sb, (m1 + m2)^2, Infinity}, 
                        Method -> {"GaussKronrodRule", "SymbolicProcessing" -> 0},
                        MaxRecursion -> 20, 
                        AccuracyGoal -> 5];
 
RNIntBBviaDudek[s_][m1_, m2_, mRecoil_, f2R_, m2R_, s0_, orbL_, scalingFac_] :=
  scalingFac /\[Pi] NIntegrate[ f2R ImagSigmads2Dudek[sb,m1,m2,f2R, m2R, s0, orbL] cmReid[s,Sqrt[sb], mRecoil]
                        , {sb, (m1 + m2)^2, Infinity}, 
                        Method -> {"GaussKronrodRule", "SymbolicProcessing" -> 0},
                        MaxRecursion -> 20, 
                        AccuracyGoal -> 5];

RNIntBBviaBW[s_][m1_, m2_, mrec_, orbMom_, M0_, Width_, scale_] :=
  scale /\[Pi] NIntegrate[ rho[sb, m1, m2] Abs[BreitWignerBlattWRelNorm[sb, orbMom, M0, Width, m1, m2]]^2 cmReid[s, Sqrt[sb], mrec]
                        , {sb, (m1 + m2)^2, Infinity}, 
                        Method -> {"GaussKronrodRule", "SymbolicProcessing" -> 0},
                        MaxRecursion -> 20, 
                        AccuracyGoal -> 5];




RNIntDudekDouble[m11_, m12_, f2R1_, m2R1_, orbL1_, s01_, m21_, m22_, f2R2_, m2R2_, orbL2_, s02_, scalingFac_, smax_] :=
  scalingFac /\[Pi] NIntegrate[ f2R1 ImagSigmads2Dudek[sb1,m11,m12,f2R1,m2R1, s01, orbL1]
				f2R2 ImagSigmads2Dudek[sb2,m21,m22,f2R2,m2R2, s02, orbL2] 
				, {sb1, (m11 + m12)^2, smax}, {sb2, (m21 + m22)^2, smax}, 
			Method -> {"GaussKronrodRule", "SymbolicProcessing" -> 0},
			MaxRecursion -> 20, 
			AccuracyGoal -> 5];

RNIntDudekDoubleViaDudek[s_][m11_, m12_, f2R1_, m2R1_, orbL1_, s01_, m21_, m22_, f2R2_, m2R2_, orbL2_, s02_, scalingFac_] :=
  scalingFac /\[Pi] NIntegrate[ f2R1 ImagSigmads2Dudek[sb1,m11,m12,f2R1,m2R1, s01, orbL1]
				f2R2 ImagSigmads2Dudek[sb2,m21,m22,f2R2,m2R2, s02, orbL2]
				cmDudek[s, Sqrt[sb1], Sqrt[sb2]] 
				, {sb1, (m11 + m12)^2, Infinity}, {sb2, (m21 + m22)^2, Infinity}, 
			Method -> {"GaussKronrodRule", "SymbolicProcessing" -> 0},
			MaxRecursion -> 20, 
			AccuracyGoal -> 5];




RNIntBBDouble[m11_, m12_, f2R1_, m2R1_, orbL1_, s01_, m21_, m22_, f2R2_, m2R2_, orbL2_, s02_, scalingFac_, smax_] :=
  scalingFac /\[Pi] NIntegrate[ f2R1 ImagSigmads2BB[sb1,m11,m12,f2R1,m2R1, s01, orbL1]
				f2R2 ImagSigmads2BB[sb2,m21,m22,f2R2,m2R2, s02, orbL2] 
				, {sb1, (m11 + m12)^2, smax}, {sb2, (m21 + m22)^2, smax}, 
			Method -> {"GaussKronrodRule", "SymbolicProcessing" -> 0},
			MaxRecursion -> 20, 
			AccuracyGoal -> 5];


RNIntDudekDoubleViaBB[s_][m11_, m12_, f2R1_, m2R1_, orbL1_, s01_, m21_, m22_, f2R2_, m2R2_, orbL2_, s02_, scalingFac_] :=
  scalingFac /\[Pi] NIntegrate[ f2R1 ImagSigmads2BB[sb1,m11,m12,f2R1,m2R1, s01, orbL1]
				f2R2 ImagSigmads2BB[sb2,m21,m22,f2R2,m2R2, s02, orbL2]
				cmDudek[s, Sqrt[sb1], Sqrt[sb2]] 
				, {sb1, (m11 + m12)^2, Infinity}, {sb2, (m21 + m22)^2, Infinity}, 
			Method -> {"GaussKronrodRule", "SymbolicProcessing" -> 0},
			MaxRecursion -> 20, 
			AccuracyGoal -> 5];

RNIntBWDouble[M01_, Width1_, m11_, m12_, orbMom1_, M02_, Width2_, m21_, m22_, orbMom2_,scalingFac_, smax_] :=
  scalingFac /\[Pi] NIntegrate[ rho[sb1, m11, m12] Abs[BreitWignerBlattWRelNorm[sb1, orbMom1, M01, Width1, m11, m12]]^2
				rho[sb2, m21, m22] Abs[BreitWignerBlattWRelNorm[sb2, orbMom2, M02, Width2, m21, m22]]^2
				, {sb1, (m11 + m12)^2, smax}, {sb2, (m21 + m22)^2, smax}, 
				Method -> {"GaussKronrodRule", "SymbolicProcessing" -> 0},
				MaxRecursion -> 20, 
				AccuracyGoal -> 5];


RNIntDudekDoubleViaBW[s_][M01_, Width1_, m11_, m12_, orbMom1_, M02_, Width2_, m21_, m22_, orbMom2_,scalingFac_] :=
  scalingFac /\[Pi] NIntegrate[ rho[sb1, m11, m12] Abs[BreitWignerBlattWRelNorm[sb1, orbMom1, M01, Width1, m11, m12]]^2
				rho[sb2, m21, m22] Abs[BreitWignerBlattWRelNorm[sb2, orbMom2, M02, Width2, m21, m22]]^2
				cmDudek[s, Sqrt[sb1], Sqrt[sb2]] 
				, {sb1, (m11 + m12)^2, 10000.}, {sb2, (m21 + m22)^2, 10000.}, 
				Method -> {"GaussKronrodRule", "SymbolicProcessing" -> 0},
				MaxRecursion -> 20, 
				AccuracyGoal -> 5]; 
