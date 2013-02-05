// epemReaction class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/DataUtils.hh"
#include "Utils/PawianCollectionUtils.hh"

class Particle;
class IsobarLSDecay;
class IsobarHeliDecay;

class epemReaction {

public:
  epemReaction(std::vector<std::pair<Particle*, Particle*> >& prodPairs);

  virtual ~epemReaction();

  virtual void print(std::ostream& os) const;
  std::vector< boost::shared_ptr<IsobarLSDecay> >& productionCanoDecays() {return _prodCanoDecs;}
  std::vector< boost::shared_ptr<IsobarHeliDecay> >& productionHeliDecays() {return _prodHeliDecs;}  

protected:

private:
  boost::shared_ptr<const jpcRes> _epemJPC;

  std::vector< boost::shared_ptr<IsobarLSDecay> > _prodCanoDecs;
  std::vector< boost::shared_ptr<IsobarHeliDecay> > _prodHeliDecs;

};
