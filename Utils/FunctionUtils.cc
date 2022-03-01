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

  Id1StringType spin1Index(const Spin& spin1){//only for spin<=9 particles
    Id1StringType result=spin1.ToIndex();
    return result;
  }

 
  Id2StringType spin2Index(const Spin& spin1, const Spin& spin2){//only for spin<=9 particles
    Id2StringType result=spin1.ToIndex()*36+spin2.ToIndex();
    return result;
  }

  Id3StringType spin3Index(const Spin& spin1, const Spin& spin2, const Spin& spin3){ //only for spin<=9 particles
    Id3StringType result=spin1.ToIndex()*36*36+spin2.ToIndex()*36+spin3.ToIndex();
    return result;
  }

  short spin1IdIndex(unsigned short id, const Spin& spin1){//only for spin<=9 particles and ids<500
    if(id>500){
      std::cout << "id value: " << id << " is not less than 500!!!" << std::endl;
      exit(1);
    }
    short result=id*72+spin1.ToIndex();
    return result;
  }

}
