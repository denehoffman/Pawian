// OmegaTo3PiLSDecay class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <algorithm>

#include "PwaUtils/OmegaTo3PiLSDecay.hh"
#include "PwaUtils/AbsEnv.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/FunctionUtils.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/EvtDataBaseList.hh"

OmegaTo3PiLSDecay::OmegaTo3PiLSDecay(Particle* mother, Particle* daughter1, Particle* daughter2, Particle* daughter3, AbsEnv* theEnv) :
  AbsDecay(mother, daughter1, daughter2, theEnv)
  ,_daughter3(daughter3)
{
  _finalStateParticles.push_back(daughter3);
  _finalStateParticlesDaughter3.push_back(daughter3);

  pawian::Collection::PtrLess thePtrLess;
  std::sort(_finalStateParticles.begin(), _finalStateParticles.end(), thePtrLess);

  _name+="_"+daughter3->name();
  _fitParamSuffix=_name;
  
  //check correct quantum numbers
  //...
  Spin validL=1;
  Spin validS=0;
  boost::shared_ptr<const JPCLS> theValidJPCLS(new JPCLS(_motherJPCPtr, validL, validS));
  _JPCLSDecAmps.push_back(theValidJPCLS);

  _wignerDKey="normOmega_"+_motherJPCPtr->name()+FunctionUtils::particleListName(_finalStateParticles);
  _lambdaDecKey="lambdaOmega_"+_motherJPCPtr->name()+FunctionUtils::particleListName(_finalStateParticles); 
}

OmegaTo3PiLSDecay::~OmegaTo3PiLSDecay(){
}

void OmegaTo3PiLSDecay::fillWignerDs(std::map<std::string , Vector4<double> >& fsMap, EvtData* evtData){
  Vector4<double> all4Vec(0.,0.,0.,0.);
  Vector4<double> mother4Vec(0.,0.,0.,0.);

  //fill all4Vec
  std::map<std::string, Vector4<double> >::iterator itMap;
  for(itMap=fsMap.begin(); itMap!=fsMap.end(); ++itMap){
    all4Vec+=itMap->second;
  }

  //fill mother4Vec
  std::vector<Particle*>::iterator itP;
  for(itP = _finalStateParticles.begin(); itP != _finalStateParticles.end(); ++itP){
    itMap=fsMap.find((*itP)->name());
    mother4Vec+=itMap->second;
  }

  //fill daughter1 and daughter2 4Vec
  itMap=fsMap.find(_daughter1->name());
  Vector4<double> daughter1_4Vec=itMap->second;
  itMap=fsMap.find(_daughter2->name());
  Vector4<double> daughter2_4Vec=itMap->second;
  itMap=fsMap.find(_daughter3->name());
  Vector4<double> daughter3_4Vec=itMap->second;

  Vector4<double> daughter1_HeliOmega=helicityVec(all4Vec, mother4Vec, daughter1_4Vec);
  Vector4<double> daughter2_HeliOmega=helicityVec(all4Vec, mother4Vec, daughter2_4Vec);
  Vector4<double> daughter3_HeliOmega=helicityVec(all4Vec, mother4Vec, daughter3_4Vec);

  //calculate normal of the decay plane
  Vector4<float> normOmegaDecHeliOmega_4V(0.5*(daughter1_HeliOmega.T()+daughter2_HeliOmega.T()+daughter3_HeliOmega.T()),
					  daughter1_HeliOmega.Y()*daughter2_HeliOmega.Z()-daughter1_HeliOmega.Z()*daughter2_HeliOmega.Y(),
					  daughter1_HeliOmega.Z()*daughter2_HeliOmega.X()-daughter1_HeliOmega.X()*daughter2_HeliOmega.Z(),        
					  daughter1_HeliOmega.X()*daughter2_HeliOmega.Y()-daughter1_HeliOmega.Y()*daughter2_HeliOmega.X());

  double theQ=daughter1_HeliOmega.E()-daughter1_HeliOmega.M()+daughter2_HeliOmega.E()-daughter2_HeliOmega.M()+daughter3_HeliOmega.E()-daughter3_HeliOmega.M();
  double lambdaNorm=theQ*theQ*(theQ*theQ/108.+daughter1_HeliOmega.M()*theQ/9.+daughter1_HeliOmega.M()*daughter1_HeliOmega.M()/3.);
  double lambdaOmegaDec=normOmegaDecHeliOmega_4V.P()*normOmegaDecHeliOmega_4V.P()/lambdaNorm;
  evtData->DoubleString["lamOmegaDec"]=lambdaOmegaDec;

  for(Spin lamOmega=-1;  lamOmega<=1; ++lamOmega){
    evtData->WignerDsString[_wignerDKey][1][lamOmega][0]=Wigner_D(normOmegaDecHeliOmega_4V.Phi(), normOmegaDecHeliOmega_4V.Theta(),0, 1,lamOmega,0);
  }
}

void OmegaTo3PiLSDecay::print(std::ostream& os) const{
  os << "\nJPCLS amplitudes for decay\t" << _name << ":\n";
  os << "suffix for fit parameter name:\t" << _fitParamSuffix << "\n";
  
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator it;
  for (it = _JPCLSDecAmps.begin(); it!= _JPCLSDecAmps.end(); ++it){
    (*it)->print(os);
    os << "\n";
  }

  AbsDecay::print(os);  
  os << "\n";
}
