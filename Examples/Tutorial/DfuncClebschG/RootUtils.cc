//************************************************************************//
//                                                                        //
//  Copyright 2016 Bertram Kopf (bertram@ep1.rub.de)                      //
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

// RootUtils class definition file. -*- C++ -*-
// Copyright 2016 Bertram Kopf

#include "Examples/Tutorial/DfuncClebschG/RootUtils.hh"
#include "TVector3.h"
#include "Utils/PawianConstants.hh"

TLorentzVector RootUtils::heliVec(const TLorentzVector& motherRef,  const TLorentzVector& ref,const TLorentzVector& mother, const TLorentzVector& daughter){
  TLorentzVector result=daughter;
  TLorentzVector refTrafo=ref;
  TLorentzVector refRecoilTrafo=motherRef-ref;

  //boost all vectors into the mother rest frame
  result.Boost(-mother.BoostVector());
  refTrafo.Boost(-mother.BoostVector());
  refRecoilTrafo.Boost(-mother.BoostVector());

  //rotate vectors so that refRecoilTrafo moves in the negative direction of the z-axis   
  result.RotateZ(-refTrafo.Phi());
  result.RotateY(PawianConstants::pi-refTrafo.Theta());

  refRecoilTrafo.RotateZ(-refTrafo.Phi());
  refRecoilTrafo.RotateY(PawianConstants::pi-refTrafo.Theta());

  //rotate around the z-axis so that refRecoil lies in the x-z plain
  result.RotateZ(PawianConstants::pi-refRecoilTrafo.Phi()); 

  return result;
}

