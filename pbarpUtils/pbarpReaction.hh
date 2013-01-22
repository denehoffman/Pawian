// pbarpReaction class definition file. -*- C++ -*-
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

class pbarpReaction {

public:
  pbarpReaction(std::vector<std::pair<Particle*, Particle*> >& prodPairs, int lmax=3);

  virtual ~pbarpReaction();

  virtual void print(std::ostream& os) const;
  std::vector< boost::shared_ptr<IsobarLSDecay> >& productionDecays() {return _prodDecs;}  
  std::vector< boost::shared_ptr<const jpcRes> >& jpcStates() {return _pbarpJPCs;}
  std::vector< boost::shared_ptr<const JPCLS> >& jpclsStates() {return _pbarpJPCLSs;}
  std::vector< boost::shared_ptr<const JPCLS> >& jpclsSingletStates() {return _pbarpJPCLSsinglet;}
  std::vector< boost::shared_ptr<const JPCLS> >& jpclsTriplet0States() {return _pbarpJPCLStriplet0;}
  std::vector< boost::shared_ptr<const JPCLS> >& jpclsTripletp1States() {return _pbarpJPCLStripletp1;}
  std::vector< boost::shared_ptr<const JPCLS> >& jpclsTripletm1States() {return _pbarpJPCLStripletm1;}

  std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<IsobarLSDecay> >, pawian::Collection::SharedPtrLess >& singletDecMap() {return _pbarpSingletDecMap;}
protected:

private:
  int _lmax;
  std::vector< boost::shared_ptr<const jpcRes> > _pbarpJPCs;
  std::vector< boost::shared_ptr<const JPCLS> > _pbarpJPCLSs;
  std::vector< boost::shared_ptr<const JPCLS> > _pbarpJPCLSsinglet;
  std::vector< boost::shared_ptr<const JPCLS> > _pbarpJPCLStriplet0;
  std::vector< boost::shared_ptr<const JPCLS> > _pbarpJPCLStripletp1;
  std::vector< boost::shared_ptr<const JPCLS> > _pbarpJPCLStripletm1;

  std::vector< boost::shared_ptr<IsobarLSDecay> > _prodDecs;

  std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<IsobarLSDecay> >, pawian::Collection::SharedPtrLess > _pbarpSingletDecMap;
  std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<IsobarLSDecay> >, pawian::Collection::SharedPtrLess > _pbarpTriplet0DecMap;
  std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<IsobarLSDecay> >, pawian::Collection::SharedPtrLess > _pbarpTripletp1DecMap;
  std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<IsobarLSDecay> >, pawian::Collection::SharedPtrLess > _pbarpTripletm1DecMap;

  void  fillMap(std::vector< boost::shared_ptr<const JPCLS> >& pbarpLSs, std::vector<boost::shared_ptr<IsobarLSDecay> >& decs, std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<IsobarLSDecay> >, pawian::Collection::SharedPtrLess > toFill);
};
