// pbarpCanoLh class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "pbarpUtils/pbarpCanoLh.hh"
#include "pbarpUtils/pbarpEnv.hh"
#include "pbarpUtils/pbarpReaction.hh"
#include "PwaUtils/LSDecAmps.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


pbarpCanoLh::pbarpCanoLh(boost::shared_ptr<const EvtDataBaseList> theEvtList) :
  pbarpBaseLh(theEvtList)
{
  initialize();
}



pbarpCanoLh::~pbarpCanoLh()
{;
}


void pbarpCanoLh::print(std::ostream& os) const{

}


void  pbarpCanoLh::initialize(){

  std::vector< boost::shared_ptr<IsobarLSDecay> > theDecs = _pbarpReactionPtr->productionDecays();
  std::vector< boost::shared_ptr<IsobarLSDecay> >::iterator it;
  for (it=theDecs.begin(); it!=theDecs.end(); ++it){
    boost::shared_ptr<AbsDecay> currentDec( (*it).get() );
    boost::shared_ptr<AbsXdecAmp> currentAmp=XdecAmpRegistry::instance()->getXdecAmp(currentDec);
    _decAmps.push_back(currentAmp);
  }

  std::vector< boost::shared_ptr<const JPCLS> > jpclsSingletStates=_pbarpReactionPtr->jpclsSingletStates();
  fillMap(jpclsSingletStates, _decAmps, _decAmpsSinglet);

  std::vector< boost::shared_ptr<const JPCLS> > jpclsTriplet0States=_pbarpReactionPtr->jpclsTriplet0States();
  fillMap(jpclsTriplet0States, _decAmps, _decAmpsTriplet0);

  std::vector< boost::shared_ptr<const JPCLS> > jpclsTripletp1States=_pbarpReactionPtr->jpclsTripletp1States();
  fillMap(jpclsTripletp1States, _decAmps, _decAmpsTripletp1);

  std::vector< boost::shared_ptr<const JPCLS> > jpclsTripletm1States=_pbarpReactionPtr->jpclsTripletm1States();
  fillMap(jpclsTripletm1States, _decAmps, _decAmpsTripletm1);  
}


