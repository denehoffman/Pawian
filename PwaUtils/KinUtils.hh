#ifndef KinUtils_HH
#define KinUtils_HH

#include <cassert>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

template<typename T>

Vector4<T> helicityVec(const Vector4<T>& cm, const Vector4<T>& mother, const Vector4<T>& daughter){
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

#endif
