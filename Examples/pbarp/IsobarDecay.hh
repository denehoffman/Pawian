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
//#include "Particle/Particle.hh"

class Particle;

class IsobarDecay {

public:
  IsobarDecay(Particle* mother, Particle* daughter1, Particle* daughter2);
  IsobarDecay(boost::shared_ptr<const jpcRes> motherJPCPtr, Particle* daughter1, Particle* daughter2, std::string motherName="pbarp");
  virtual ~IsobarDecay();
  //  virtual IsobarDecay* clone_() const = 0;
  virtual const std::string name() const {return _name;} 
  virtual std::string fitParSuffix() const {return _fitParamSuffix;}
  void setFitParSuffix(std::string& suffix) {_fitParamSuffix = suffix;}
  boost::shared_ptr<const jpcRes> motherJPC(){ return _motherJPCPtr;}
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSAmps(){ return _JPCLSDecAmps;}
  virtual void print(std::ostream& os) const;
  
  Particle* motherPart() {return _mother;}

protected:
  Particle* _mother;
  Particle* _daughter1;
  Particle* _daughter2;

  boost::shared_ptr<const jpcRes> _motherJPCPtr;
  boost::shared_ptr<const jpcRes> _daughter1JPCPtr;
  boost::shared_ptr<const jpcRes> _daughter2JPCPtr;

  std::string _name;
  std::string _fitParamSuffix;

  std::vector< boost::shared_ptr<const JPCLS> > _JPCLSDecAmps;

  boost::shared_ptr<IsobarDecay> _isoDecDaughter1;
  boost::shared_ptr<IsobarDecay> _isoDecDaughter2;  
};
