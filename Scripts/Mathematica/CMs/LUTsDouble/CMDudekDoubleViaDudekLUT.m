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

Get[FileNameJoin[{Directory[], "../Integrals.m"}]]
(* Get[FileNameJoin[{Directory[], "../InputParamsDoubleSigmaSigma.m"}]] *)
Get[FileNameJoin[{Directory[], "../InputParamsDoubleOmegaPhi.m"}]]

(* Precision[1.2]; *)
(* Precision[12/10]; *)



mrealstepsize:= 0.0005;
mrealmin:= 0.026;
mrealmax:= 2.8;

mimagstepsize := 0.01;
numOfCores := 8;
mimagstepsizeNoCtimes := numOfCores mimagstepsize;
mimagmin := -0.5;
mimagmax := 0.01;
 
(*
{t1,b1} =
 AbsoluteTiming[
		*)		
For[ mreal=mrealmin, mreal<=mrealmax, mreal+=mrealstepsize,
     Print["mreal: ", mreal];
     ofile=OpenAppend["./CMDudekDoubleViaDudekLUT.dat"];
     For[ mimag = mimagmax, mimag >= mimagmin, mimag -= mimagstepsizeNoCtimes,{
	  SetSharedVariable[list];
          list = {};
	  ParallelDo[{
	      currentCM = RNIntDudekDoubleViaDudek[(mreal + I mimag) (mreal + I mimag)][mDec11, mDec12, f2PoleDudek1, m2PoleDudek1, orbMomBW1, (mDec11+mDec12)^2, mDec21, mDec22, f2PoleDudek2, m2PoleDudek2, orbMomBW2,(mDec21+mDec22)^2 , scalingFacDudekSmearing];
	      AppendTo[list, {mimag, currentCM}]},
	    {mimag, mimag, mimag-mimagstepsizeNoCtimes+mimagstepsize, -mimagstepsize}
	    ],
	  list = SortBy[list, First, Greater];
	  WriteString[ofile, NumberForm[mreal,5], "\t", NumberForm[#[[1]], 5], "\t", NumberForm[Re[#[[2]]], 16], "\t", NumberForm[Im[#[[2]]], 16], "\n"] & /@ list;
       };
       ]
     Close[ofile];
     ];
 (* 
];

Print["t1: ", t1,"\tb1: ", b1];
 *)
