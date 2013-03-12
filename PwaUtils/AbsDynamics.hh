// AbsDynamics class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <mutex>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "PwaUtils/AbsParamHandler.hh"

class Particle;
class AbsDecay;
class AbsXdecAmp;

class AbsDynamics : public AbsParamHandler{

public:
  AbsDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother );
  virtual ~AbsDynamics();

  virtual std::string name() {return _name;};
  virtual complex<double> eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom=0)=0;
  virtual void cacheAmplitudes();
  virtual void fillMasses(EvtData* theData);
  virtual void setMassKey(std::string& theMassKey){_massKey=theMassKey;}
  virtual void addGrandMa(boost::shared_ptr<AbsDecay> theDec) {;}
  virtual const std::string& grandMaKey(AbsXdecAmp* grandmaAmp){return _grandmaKey;}

protected:
  std::string _name;
  std::string _massKey;
  std::vector<Particle*> _fsParticles;
  Particle* _mother;
  std::map<int, complex<float> >  _cachedMap;
  std::mutex theMutex;
  std::string _dynKey;
  std::string _grandmaKey;
private:

};
