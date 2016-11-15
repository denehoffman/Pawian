//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
//									  //
//  This file is part of Pawian.					  //
//									  //
//  Pawian is free software: you can redistribute it and/or modify	  //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	  //
//									  //
//  Pawian is distributed in the hope that it will be useful,		  //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  //
//  GNU General Public License for more details.	      		  //
//									  //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	  //
//									  //
//************************************************************************//

#pragma once

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

namespace KinUtils{
  template<typename T> 
  Vector4<T> heliVec(const Vector4<T>& motherRef,  const Vector4<T>& ref,const Vector4<T>& mother, const Vector4<T>& daughter);
  
  template<typename T>   
  Vector4<T> gottfriedJacksonVec(const Vector4<T>& motherRef,  const Vector4<T>& ref,const Vector4<T>& mother, const Vector4<T>& daughter);
  
  template<typename T>
  Vector4<T> perpTo(const Vector4<T>& vecA,  const Vector4<T>& vecB);

  void SetTheta(Vector4<double>& the4Vec, double theta);
  void SetPhi(Vector4<double>& the4Vec, double phi);
  Vector4<double> heliEvtGenVec(Vector4<double>& motherRef, Vector4<double>& ref, Vector4<double>& mother, Vector4<double>& daughter); 
}


template<typename T>
Vector4<T> helicityVec(const Vector4<T>& cm, const Vector4<T>& mother, 
		       const Vector4<T>& daughter) {
  Vector4<T> mother_cm(mother);
  mother_cm.Boost(cm);
  
  Vector4<T> mother_cm_rot(mother_cm);
  mother_cm_rot.RotateZ(-mother_cm.Phi());
  mother_cm_rot.RotateY(-mother_cm.Theta()); 

  Vector4<T> result(daughter);
  result.Boost(cm);
  result.RotateZ(-mother_cm.Phi());
  result.RotateY(-mother_cm.Theta()); 
  
  result.Boost(mother_cm_rot);

  return result;
}

template<typename T>
T costDecHeli(const Vector4<T>& parent,const Vector4<T>& resonance,
	      const Vector4<T>& daughter) {
  T pd=parent*daughter;
  T pr=parent*resonance;
  T rd=resonance*daughter;
  T mp2=parent.Mass2();
  T mr2=resonance.Mass2();
  T md2=daughter.Mass2();

  T cost=(pd*mr2-pr*rd)/sqrt((pr*pr-mr2*mp2)*(rd*rd-mr2*md2));

  return cost;
}
