#include <algorithm>   

#include "Utils/FunctionUtils.hh"
#include "Particle/Particle.hh"
#include "qft++/relativistic-quantum-mechanics/Spin.hh"

namespace FunctionUtils{


  std::string particleListName(std::vector<Particle*>& particleVec){
  std::string result;

  std::vector<std::string> nameVec;
  std::vector<Particle*>::iterator it;
  for(it=particleVec.begin(); it!=particleVec.end(); ++it){
    nameVec.push_back((*it)->name());
  }

  std::sort(nameVec.begin(), nameVec.end());
  std::vector<std::string>::iterator itStr;
  for(itStr=nameVec.begin(); itStr!=nameVec.end(); ++itStr){
    result+=(*itStr);
  }

  return result;
  }
 
  Id2StringType spin2Index(Spin& spin1, Spin& spin2){//only for spin<=9 particles
    Id2StringType result=spin1.ToIndex()*36+spin2.ToIndex();
    return result;
  }

  Id3StringType spin3Index(Spin& spin1, Spin& spin2, Spin& spin3){ //only for spin<=9 particles
    Id3StringType result=spin1.ToIndex()*36*36+spin2.ToIndex()*36+spin3.ToIndex();
    return result;
  }
}
