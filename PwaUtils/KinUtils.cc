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

// DataUtils class definition file. -*- C++ -*-
// Copyright 2016 Bertram Kopf

#include "PwaUtils/KinUtils.hh"
#include "PspGen/EvtVector4R.hh"
#include "PspGen/EvtKine.hh"
#include "Utils/PawianConstants.hh"

void KinUtils::SetTheta(Vector4<double>& the4Vec, double theta){
  double energy = the4Vec.E();
  double mom = the4Vec.P();
  double phi   = the4Vec.Phi();
  
  double px = mom*sin(theta)*cos(phi);
  double py = mom*sin(theta)*sin(phi);
  double pz = mom*cos(theta);
  
  the4Vec.SetV4(energy, px, py, pz);
}

void KinUtils::SetPhi(Vector4<double>& the4Vec, double phi) 
{
  double p_xy   = the4Vec.Pxy();
  
  double energy = the4Vec.E();
  double px = p_xy*cos(phi);
  double py = p_xy*sin(phi);      
  double pz = the4Vec.Pz();
  
  the4Vec.SetV4(energy, px, py, pz);
}

Vector4<double>  KinUtils::heliEvtGenVec(Vector4<double>& motherRef, Vector4<double>& ref, Vector4<double>& mother, 
				   Vector4<double>& daughter) {
  //use EvtGen classes for calculating theta and phi; temporary solution and should be fixed soon
  EvtVector4R motherRefEvtGen(motherRef.E(), motherRef.Px(), motherRef.Py(), motherRef.Pz());  
  EvtVector4R refEvtGen(ref.E(), ref.Px(), ref.Py(), ref.Pz());
  EvtVector4R motherEvtGen(mother.E(), mother.Px(), mother.Py(), mother.Pz());
  EvtVector4R daughterEvtGen(daughter.E(), daughter.Px(), daughter.Py(), daughter.Pz()); 
  
  double thetaHeliDaughter= acos(EvtDecayAngle(refEvtGen, motherEvtGen, daughterEvtGen));
  double phiHeliDaughter= EvtDecayAnglePhi(motherRefEvtGen, refEvtGen, motherEvtGen, daughterEvtGen);
  
  Vector4<double> result(daughter);
  result.Boost(mother);
  KinUtils::SetTheta(result, thetaHeliDaughter);
  KinUtils::SetPhi(result, phiHeliDaughter);
  
  return result;
}

template<typename T> 
Vector4<T> KinUtils::heliVec(const Vector4<T>& motherRef,  const Vector4<T>& ref,const Vector4<T>& mother, const Vector4<T>& daughter){
  Vector4<T> result=daughter;
  Vector4<T> refTrafo=ref;
  Vector4<T> refRecoilTrafo=motherRef-ref;

  //boost all vectors into the mother rest frame
  result.Boost(mother);
  refTrafo.Boost(mother);
  refRecoilTrafo.Boost(mother);

  //rotate vectors so that refRecoilTrafo moves in the negative direction of the z-axis   
  result.RotateZ(-refTrafo.Phi());
  result.RotateY(PawianConstants::pi-refTrafo.Theta());

  refRecoilTrafo.RotateZ(-refTrafo.Phi());
  refRecoilTrafo.RotateY(PawianConstants::pi-refTrafo.Theta());

  //rotate around the z-axis so that refRecoil lies in the x-z plain
  result.RotateZ(PawianConstants::pi-refRecoilTrafo.Phi()); 
  
  return result;
}

template<typename T>
Vector4<T> KinUtils::perpTo(const Vector4<T>& vecA,  const Vector4<T>& vecB){

  Vector4<T> result(vecA.E()+vecB.E(),
		    vecA.Y()*vecB.Z()-vecA.Z()*vecB.Y(),
		    vecA.Z()*vecB.X()-vecA.X()*vecB.Z(),
		    vecA.X()*vecB.Y()-vecA.Y()*vecB.X());
  // energy is arbitrary; use sum of vector A and vector B here
 return result;
} 


template Vector4<double> KinUtils::heliVec(const Vector4<double>&,  const Vector4<double>&,const Vector4<double>& , const Vector4<double>& );  
template Vector4<float> KinUtils::heliVec(const Vector4<float>&,  const Vector4<float>&,const Vector4<float>& , const Vector4<float>& );
template Vector4<double> KinUtils::perpTo(const Vector4<double>& vecA,  const Vector4<double>& vecB);
template Vector4<float> KinUtils::perpTo(const Vector4<float>& vecA,  const Vector4<float>& vecB);
