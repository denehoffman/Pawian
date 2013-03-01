// epemBaseLh class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "epemUtils/epemBaseLh.hh"
#include "epemUtils/epemEnv.hh"
#include "epemUtils/epemReaction.hh"
#include "PwaUtils/HeliDecAmps.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


epemBaseLh::epemBaseLh(boost::shared_ptr<const EvtDataBaseList> theEvtList) :
  AbsLh(theEvtList)
  ,_usePhasespace(false)
  ,_phasespaceKey("Phasespace")
  ,_highestJFsp(0)
  ,_isHighestJaPhoton(true)
{
  initialize();
}



epemBaseLh::~epemBaseLh()
{;
}


complex<double> epemBaseLh::calcSpinDensity(Spin M1, Spin M2, std::string& nameDec, EvtData* theData){

  complex<double> result(0.,0.);

  return result;
}




complex<double> epemBaseLh::calcProdPartAmp(Spin lamX, Spin lamDec, std::string nameDec, EvtData* theData, 
					     std::map <boost::shared_ptr<const JPCLS>,
					     std::vector< boost::shared_ptr<AbsXdecAmp> >,
					     pawian::Collection::SharedPtrLess > pbarpAmps){
   complex<double> resultAmp(0.,0.);

   return resultAmp;
}


double epemBaseLh::calcEvtIntensity(EvtData* theData, fitParams& theParamVal){

  double result=0.;

  int lamSteps=1;
  if(_isHighestJaPhoton) lamSteps=2;

  for (int lamHigestJFsp=-_highestJFsp; lamHigestJFsp<=_highestJFsp; lamHigestJFsp=lamHigestJFsp+lamSteps){
 
  complex<double> lamp1Amp(0.,0.);

  std::vector<boost::shared_ptr<AbsXdecAmp> >::iterator itDec;

  for( itDec=_decAmps.begin(); itDec!=_decAmps.end(); ++itDec){
    complex<double> currentDecAmp=(*itDec)->XdecAmp(1, theData, lamHigestJFsp);
    lamp1Amp+=currentDecAmp;
  }
  

 complex<double> lamm1Amp(0.,0.);
 for( itDec=_decAmps.begin(); itDec!=_decAmps.end(); ++itDec){
   complex<double> currentDecAmp=(*itDec)->XdecAmp(-1, theData, lamHigestJFsp);
   lamm1Amp+=currentDecAmp;
 }

  result += norm(lamp1Amp) + norm(lamm1Amp); 

  }


  if(_usePhasespace){
    result = result + theParamVal.otherParams[_phasespaceKey];
  }
  return result;  

}

void epemBaseLh::print(std::ostream& os) const{

}


void  epemBaseLh::initialize(){

  std::vector<Particle*> fsParticles=epemEnv::instance()->finalStateParticles();
  std::vector<Particle*>::iterator itParticle;
  bool highestJFound=false;

  for (itParticle=fsParticles.begin(); itParticle != fsParticles.end(); ++itParticle){
    int current2J = (*itParticle)->twoJ();
    if(current2J>0){
      if(highestJFound){
	Alert << "final states with more than 1 particles with J>0 not supported!!!!" << endmsg;
	exit(1); 
      }
      _highestJFsp=int(current2J/2);
      if( (*itParticle)->name() != "photon" ) _isHighestJaPhoton=false;
    }
  }  
  
  _epemReactionPtr = epemEnv::instance()->reaction();

std::vector< boost::shared_ptr<IsobarHeliDecay> > theDecs = _epemReactionPtr->productionHeliDecays();
  std::vector< boost::shared_ptr<IsobarHeliDecay> >::iterator it;
  for (it=theDecs.begin(); it!=theDecs.end(); ++it){
    //    boost::shared_ptr<AbsDecay> currentDec((*it).get() );
    //    boost::shared_ptr<AbsXdecAmp> currentAmp=XdecAmpRegistry::instance()->getXdecAmp(currentDec);
    boost::shared_ptr<AbsXdecAmp> currentAmp=XdecAmpRegistry::instance()->getXdecAmp((*it)->absDecPtr());
    _decAmps.push_back(currentAmp);
  }

}




