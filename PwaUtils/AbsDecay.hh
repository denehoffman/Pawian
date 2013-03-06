// AbsDecay class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#include "PwaUtils/DataUtils.hh"
#include "Utils/PawianCollectionUtils.hh"

class Particle;
class EvtData;
class AbsEnv;
class AbsDynamics;

class AbsDecay : public boost::enable_shared_from_this<AbsDecay>{

public:
  AbsDecay(Particle* mother, Particle* daughter1, Particle* daughter2, AbsEnv* theEnv);
  AbsDecay(boost::shared_ptr<const jpcRes> motherJPCPtr, Particle* daughter1, Particle* daughter2, AbsEnv* theEnv, std::string motherName);
  virtual ~AbsDecay();
  //  virtual AbsDecay* clone_() const = 0;
  virtual const std::string name() const {return _name;}
  const std::string wignerDKey() {return _wignerDKey;} 
  virtual std::string fitParSuffix() const {return _fitParamSuffix;}
  void setFitParSuffix(std::string& suffix) {_fitParamSuffix = suffix;}
  virtual std::string& massParKey() {return _massParamKey;}
  void setMassParKey(const std::string& newKey) {_massParamKey = newKey;}
  boost::shared_ptr<const jpcRes> motherJPC(){ return _motherJPCPtr;}
  //  std::vector< boost::shared_ptr<const JPCLS> > JPCLSAmps(){ return _JPCLSDecAmps;}
  boost::shared_ptr<AbsDecay> decDaughter1() {return _absDecDaughter1;}
  boost::shared_ptr<AbsDecay> decDaughter2() {return _absDecDaughter2;}
  bool hasMother() {return _hasMotherPart;}
  bool isDaughter1Stable() {return _daughter1IsStable;}
  bool isDaughter2Stable() {return _daughter2IsStable;}

  std::vector<Particle*> finalStateParticles() {return _finalStateParticles;}
  std::vector<Particle*> finalStateParticlesDaughter2() {return _finalStateParticlesDaughter2;}
  virtual void fillWignerDs(std::map<std::string , Vector4<double> >& fsMap, EvtData* evtData);
  void enableDynamics(std::string& dynString, std::vector<std::string>& additionalStringVec);
  boost::shared_ptr<AbsDynamics> getDynamics(){return _absDynPtr;}
  virtual void print(std::ostream& os) const;
  
  Particle* motherPart() {return _mother;}
  Particle* daughter1Part() {return _daughter1;}
  Particle* daughter2Part() {return _daughter2;}

  virtual std::string dynType() {return _dynType;}
  virtual std::string& dynKey() {return _dynKey;}
  virtual std::pair<Particle*, Particle*>& firstDecayChannel() { return _decPair1stChannel;}
  virtual std::pair<Particle*, Particle*>& secondDecayChannel() { return _decPair2ndChannel;}
  virtual boost::shared_ptr<AbsDecay> absDecPtr() {return shared_from_this();}
  virtual std::string type() =0;

protected:
  Particle* _mother;
  Particle* _daughter1;
  Particle* _daughter2;

  bool _daughter1IsStable;
  bool _daughter2IsStable;
  bool _hasMotherPart;

  boost::shared_ptr<const jpcRes> _motherJPCPtr;
  boost::shared_ptr<const jpcRes> _daughter1JPCPtr;
  boost::shared_ptr<const jpcRes> _daughter2JPCPtr;

  std::string _name;
  std::string _fitParamSuffix;
  std::string _massParamKey;
  std::string _wignerDKey;
  std::string _dynType;
  std::string _dynKey;

  std::vector< boost::shared_ptr<const JPCLS> > _JPCLSDecAmps;

  boost::shared_ptr<AbsDecay> _absDecDaughter1;
  boost::shared_ptr<AbsDecay> _absDecDaughter2;

  std::vector<Particle*> _finalStateParticles;
  std::vector<Particle*> _finalStateParticlesDaughter1;
  std::vector<Particle*> _finalStateParticlesDaughter2;
  std::pair<Particle*, Particle*> _decPair1stChannel;
  std::pair<Particle*, Particle*> _decPair2ndChannel;

  AbsEnv* _env;
  boost::shared_ptr<AbsDynamics> _absDynPtr;
};
