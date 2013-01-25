#pragma once

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

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
Vector4<T> helicityVecInverse(const Vector4<T>& cm, const Vector4<T>& mother, 
		       const Vector4<T>& helicityVec) {
  Vector4<T> mother_cm(mother);
  mother_cm.Boost(cm);
  
  Vector4<T> mother_cm_rot(mother_cm);
  mother_cm_rot.RotateZ(-mother_cm.Phi());
  mother_cm_rot.RotateY(-mother_cm.Theta()); 

  Vector4<T> result(helicityVec);
  result.Boost(-1.* mother_cm_rot);
  
  result.RotateY(mother_cm.Theta()); 
  result.RotateZ(mother_cm.Phi());

  
  result.Boost(-1.* cm);
  
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

