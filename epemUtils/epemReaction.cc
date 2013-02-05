// epemReaction class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>

#include "epemUtils/epemReaction.hh"
#include "epemUtils/epemEnv.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"

epemReaction::epemReaction(std::vector<std::pair<Particle*, Particle*> >& prodPairs) :
  _epemJPC(new jpcRes(1,-1,-1))
{
    std::vector<std::pair<Particle*, Particle*> >::iterator itPartPairs;
    for (itPartPairs=prodPairs.begin(); itPartPairs!= prodPairs.end(); ++itPartPairs){
      //      std::string decName=(*itJPC)->name();
      boost::shared_ptr<IsobarLSDecay> currentDec(new IsobarLSDecay( _epemJPC, itPartPairs->first, itPartPairs->second, epemEnv::instance(), "epem"));

      if (currentDec->JPCLSAmps().size()>0){
	_prodCanoDecs.push_back(currentDec);
	boost::shared_ptr<IsobarHeliDecay> currentHeliDec(new IsobarHeliDecay( _epemJPC,itPartPairs->first, itPartPairs->second, epemEnv::instance(), "epem"));
	_prodHeliDecs.push_back(currentHeliDec);
      }
    }

}

epemReaction::~epemReaction(){
}

void epemReaction::print(std::ostream& os) const{
  os << "\n e+e- reaction\n";
 
  os << "\n ***** decay chains *******\n";
  std::vector< boost::shared_ptr<IsobarLSDecay> >::const_iterator itIso;
  for( itIso=_prodCanoDecs.begin(); itIso!=_prodCanoDecs.end(); ++itIso){
    (*itIso)->print(os);
  }  
}
