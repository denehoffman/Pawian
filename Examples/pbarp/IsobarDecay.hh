#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/DataUtils.hh"
#include "Utils/PawianCollectionUtils.hh"

class Particle;
class EvtDataNew;

class IsobarDecay {

public:
  IsobarDecay(Particle* mother, Particle* daughter1, Particle* daughter2);
  IsobarDecay(boost::shared_ptr<const jpcRes> motherJPCPtr, Particle* daughter1, Particle* daughter2, std::string motherName="pbarp");
  virtual ~IsobarDecay();
  //  virtual IsobarDecay* clone_() const = 0;
  virtual const std::string name() const {return _name;}
  const std::string wignerDKey() {return _wignerDKey;} 
  virtual std::string fitParSuffix() const {return _fitParamSuffix;}
  void setFitParSuffix(std::string& suffix) {_fitParamSuffix = suffix;}
  boost::shared_ptr<const jpcRes> motherJPC(){ return _motherJPCPtr;}
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSAmps(){ return _JPCLSDecAmps;}
  boost::shared_ptr<IsobarDecay> decDaughter1() {return _isoDecDaughter1;}
  boost::shared_ptr<IsobarDecay> decDaughter2() {return _isoDecDaughter2;}
  bool hasMother() {return _hasMotherPart;}
  bool isDaughter1Stable() {return _daughter1IsStable;}
  bool isDaughter2Stable() {return _daughter2IsStable;}
  bool withDynamics() {return _withDyn;}
  std::vector<Particle*> finalStateParticles() {return _finalStateParticles;}
  std::vector<Particle*> finalStateParticlesDaughter2() {return _finalStateParticlesDaughter2;}
  void fillWignerDs(std::map<std::string , Vector4<double> >& fsMap, EvtDataNew* evtData);
  void enableDynamics(std::string& dynString) {_withDyn=true;}
  virtual void print(std::ostream& os) const;
  
  Particle* motherPart() {return _mother;}
  Particle* daughter1Part() {return _daughter1;}
  Particle* daughter2Part() {return _daughter2;}

protected:
  Particle* _mother;
  Particle* _daughter1;
  Particle* _daughter2;

  bool _daughter1IsStable;
  bool _daughter2IsStable;
  bool _hasMotherPart;
  bool _withDyn;
  boost::shared_ptr<const jpcRes> _motherJPCPtr;
  boost::shared_ptr<const jpcRes> _daughter1JPCPtr;
  boost::shared_ptr<const jpcRes> _daughter2JPCPtr;

  std::string _name;
  std::string _fitParamSuffix;
  std::string _wignerDKey;

  std::vector< boost::shared_ptr<const JPCLS> > _JPCLSDecAmps;

  boost::shared_ptr<IsobarDecay> _isoDecDaughter1;
  boost::shared_ptr<IsobarDecay> _isoDecDaughter2;

  std::vector<Particle*> _finalStateParticles;
  std::vector<Particle*> _finalStateParticlesDaughter2;

};
