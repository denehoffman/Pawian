//************************************************************************//
//                                                                        //
//  Copyright 2023 Bertram Kopf (bertram@ep1.rub.de)                      //
//                 - Ruhr-Universität Bochum                              //
//                                                                        //
//  This file is part of Pawian.                                          //
//                                                                        //
//  Pawian is free software: you can redistribute it and/or modify        //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or     //
//  (at your option) any later version.                                   //
//                                                                        //
//  Pawian is distributed in the hope that it will be useful,             //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of        //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         //
//  GNU General Public License for more details.                          //
//                                                                        //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.       //
//                                                                        //
//************************************************************************//

#pragma once

#include <iostream>

const double protonMass=0.93827208816;

class KinematicsTools
{
 public:
  static double phaseSpaceFac(double mass, double massDec1, double massDec2) {
    if( massDec1<0. || massDec2<0. || mass<0.){
      std::cout << "all masses are required to be positive!!!" << std::endl;
    exit(1);
  }
  if( (massDec1+massDec2) > mass){
    std::cout << "(massDec1+massDec2) = " << massDec1+massDec2
	      << " is required to be smaller than mass = "
	      << mass << " !!!" << std::endl;
    exit(1);
  }
  double termPlus=(massDec1+massDec2)/mass;
  double termMinus=(massDec1-massDec2)/mass;
  double result=std::sqrt((1.-termPlus*termPlus) * (1.-termMinus*termMinus));
  return result;  
}
  
  static double breakupMom(double mass, double massDec1, double massDec2){
    return phaseSpaceFac(mass, massDec1, massDec2)*mass/2.;
  }

 static double pbarpMassFrompbarp(double pbarMom){
   TLorentzVector pbarpVec(0.,0.,pbarMom,sqrt(protonMass*protonMass+pbarMom*pbarMom)+protonMass);
   return pbarpVec.M();  
  }

  static double pbarMomFrompbarp(double pbarpMass){
    double pbarpSqr=std::pow((pbarpMass*pbarpMass-2.*protonMass*protonMass)/(2.*protonMass), 2.)-protonMass*protonMass;
    return sqrt( pbarpSqr);  
  }
  
};

