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
 

  unsigned int spin3Index(Spin& spin1, Spin& spin2, Spin& spin3){
    unsigned int result=spin1.ToIndex()*100*100+spin2.ToIndex()*100+spin3.ToIndex();
    return result;
  }
}
