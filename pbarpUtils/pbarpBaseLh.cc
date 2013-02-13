// pbarpBaseLh class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "pbarpUtils/pbarpBaseLh.hh"
#include "pbarpUtils/pbarpEnv.hh"
#include "pbarpUtils/pbarpReaction.hh"
#include "PwaUtils/LSDecAmps.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


pbarpBaseLh::pbarpBaseLh(boost::shared_ptr<const EvtDataBaseList> theEvtList) :
  AbsLh(theEvtList)
  ,_usePhasespace(false)
  ,_phasespaceKey("Phasespace")
  ,_highestJFsp(0)
  ,_isHighestJaPhoton(true)
{
  initialize();
}



pbarpBaseLh::~pbarpBaseLh()
{;
}


complex<double> pbarpBaseLh::calcSpinDensity(Spin M1, Spin M2, std::string& nameDec, EvtData* theData){

  complex<double> result(0.,0.);

   complex<double> singletAmpM1   = calcProdPartAmp( 0, M1, nameDec, theData, _decAmpsSinglet);
   complex<double> triplet0AmpM1  = calcProdPartAmp( 0, M1, nameDec, theData, _decAmpsTriplet0);
   complex<double> triplet1AmpM1  = calcProdPartAmp( 1, M1, nameDec, theData, _decAmpsTripletp1);
   complex<double> tripletm1AmpM1 = calcProdPartAmp(-1, M1, nameDec, theData, _decAmpsTripletm1);

   complex<double> singletAmpM2   = calcProdPartAmp( 0, M2, nameDec, theData, _decAmpsSinglet);
   complex<double> triplet0AmpM2  = calcProdPartAmp( 0, M2, nameDec, theData, _decAmpsTriplet0);
   complex<double> triplet1AmpM2  = calcProdPartAmp( 1, M2, nameDec, theData, _decAmpsTripletp1);
   complex<double> tripletm1AmpM2 = calcProdPartAmp(-1, M2, nameDec, theData, _decAmpsTripletm1);

   result = singletAmpM1   * conj(singletAmpM2) +
            triplet0AmpM1  * conj(triplet0AmpM2) +
            triplet1AmpM1  * conj(triplet1AmpM2) +
            tripletm1AmpM1 * conj(tripletm1AmpM2);

   double norm = calcSpinDensityNorm(nameDec, theData);

   return (result / norm);
}



double pbarpBaseLh::calcSpinDensityNorm(std::string& nameDec, EvtData* theData){

   double result = 0;

   for (Spin M=-1; M<=1; M++){
     complex<double> singletAmpM1   = calcProdPartAmp( 0, M, nameDec, theData, _decAmpsSinglet);
     complex<double> triplet0AmpM1  = calcProdPartAmp( 0, M, nameDec, theData, _decAmpsTriplet0);
     complex<double> triplet1AmpM1  = calcProdPartAmp( 1, M, nameDec, theData, _decAmpsTripletp1);
     complex<double> tripletm1AmpM1 = calcProdPartAmp(-1, M, nameDec, theData, _decAmpsTripletm1);

     result += (norm(singletAmpM1) + norm(triplet0AmpM1) + norm(triplet1AmpM1) + norm(tripletm1AmpM1));
   }

   return result;
}


complex<double> pbarpBaseLh::calcProdPartAmp(Spin lamX, Spin lamDec, std::string nameDec, EvtData* theData, 
					     std::map <boost::shared_ptr<const JPCLS>,
					     std::vector< boost::shared_ptr<AbsXdecAmp> >,
					     pawian::Collection::SharedPtrLess > pbarpAmps){
   complex<double> resultAmp(0.,0.);

   std::map <boost::shared_ptr<const JPCLS>,
             std::vector< boost::shared_ptr<AbsXdecAmp> >,
             pawian::Collection::SharedPtrLess >::iterator it;

   for(it=pbarpAmps.begin(); it!=pbarpAmps.end(); ++it){
     complex<double> tmpAmp(0.,0.);
     boost::shared_ptr<const JPCLS> theJPCLS=it->first;
     std::vector<boost::shared_ptr<AbsXdecAmp> > decAmps=it->second;
     std::vector<boost::shared_ptr<AbsXdecAmp> >::iterator itDec;

     for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
       Particle* particle1 = (*itDec)->absDec()->daughter1Part();
       Particle* particle2 = (*itDec)->absDec()->daughter2Part();
       complex<double> currentDecAmp;
       if(particle1->name() == nameDec)
	  currentDecAmp=(*itDec)->XdecPartAmp(lamX, lamDec, 1,  theData, 0);
       else if(particle2->name() == nameDec)
	  currentDecAmp=(*itDec)->XdecPartAmp(lamX, lamDec, 2,  theData, 0);
       else
	  Alert << "Particle " << nameDec << " not found in calcProdAmp" << endmsg;

       tmpAmp+=currentDecAmp;
     }

     double theMag=_currentParamMags[theJPCLS];
     double thePhi=_currentParamPhis[theJPCLS];
     complex<double> expi(cos(thePhi), sin(thePhi));
     tmpAmp*=theJPCLS->preFactor*theMag*expi;
     resultAmp+=tmpAmp;
   }

   return resultAmp;
}


double pbarpBaseLh::calcEvtIntensity(EvtData* theData, fitParams& theParamVal){

  double result=0.;

  std::map <boost::shared_ptr<const JPCLS>, std::vector< boost::shared_ptr<AbsXdecAmp> >, pawian::Collection::SharedPtrLess >::iterator it;

  int lamSteps=1;
  if(_isHighestJaPhoton) lamSteps=2;

  for (int lamHigestJFsp=-_highestJFsp; lamHigestJFsp<=_highestJFsp; lamHigestJFsp=lamHigestJFsp+lamSteps){
 
  complex<double> singletAmp(0.,0.);
  for(it=_decAmpsSinglet.begin(); it!=_decAmpsSinglet.end(); ++it){
    complex<double> tmpAmp(0.,0.);
    boost::shared_ptr<const JPCLS> theJPCLS=it->first;
    std::vector<boost::shared_ptr<AbsXdecAmp> > decAmps=it->second;    
     std::vector<boost::shared_ptr<AbsXdecAmp> >::iterator itDec;
     for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
       complex<double> currentDecAmp=(*itDec)->XdecAmp(0, theData, lamHigestJFsp);
       tmpAmp+=currentDecAmp;
     }

     double theMag=_currentParamMags[theJPCLS];
     double thePhi=_currentParamPhis[theJPCLS];
     complex<double> expi(cos(thePhi), sin(thePhi));
     tmpAmp*=theJPCLS->preFactor*theMag*expi;
     singletAmp+=tmpAmp;
  }

  complex<double> triplet0Amp(0.,0.);
  for(it=_decAmpsTriplet0.begin(); it!=_decAmpsTriplet0.end(); ++it){
    complex<double> tmpAmp(0.,0.);
    boost::shared_ptr<const JPCLS> theJPCLS=it->first;
    std::vector<boost::shared_ptr<AbsXdecAmp> > decAmps=it->second;    
     std::vector<boost::shared_ptr<AbsXdecAmp> >::iterator itDec;
     for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
       tmpAmp+=(*itDec)->XdecAmp(0, theData, lamHigestJFsp);
     }

     double theMag=_currentParamMags[theJPCLS];
     double thePhi=_currentParamPhis[theJPCLS];
     complex<double> expi(cos(thePhi), sin(thePhi));
     tmpAmp*=theJPCLS->preFactor*theMag*expi;
     triplet0Amp+=tmpAmp;
  }

   complex<double> tripletp1Amp(0.,0.);
  for(it=_decAmpsTripletp1.begin(); it!=_decAmpsTripletp1.end(); ++it){
    complex<double> tmpAmp(0.,0.);
    boost::shared_ptr<const JPCLS> theJPCLS=it->first;
    std::vector<boost::shared_ptr<AbsXdecAmp> > decAmps=it->second;    
     std::vector<boost::shared_ptr<AbsXdecAmp> >::iterator itDec;
     for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
       tmpAmp+=(*itDec)->XdecAmp(1, theData, lamHigestJFsp);
     }

     double theMag=_currentParamMags[theJPCLS];
     double thePhi=_currentParamPhis[theJPCLS];
     complex<double> expi(cos(thePhi), sin(thePhi));
     tmpAmp*=theJPCLS->preFactor*theMag*expi;
     tripletp1Amp+=tmpAmp;
  }  

   complex<double> tripletm1Amp(0.,0.);
  for(it=_decAmpsTripletm1.begin(); it!=_decAmpsTripletm1.end(); ++it){
    complex<double> tmpAmp(0.,0.);
    boost::shared_ptr<const JPCLS> theJPCLS=it->first;
    std::vector<boost::shared_ptr<AbsXdecAmp> > decAmps=it->second;    
     std::vector<boost::shared_ptr<AbsXdecAmp> >::iterator itDec;
     for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
       tmpAmp+=(*itDec)->XdecAmp(-1, theData, lamHigestJFsp);
     }

     double theMag=_currentParamMags[theJPCLS];
     double thePhi=_currentParamPhis[theJPCLS];
     complex<double> expi(cos(thePhi), sin(thePhi));
     tmpAmp*=theJPCLS->preFactor*theMag*expi;
     tripletm1Amp+=tmpAmp;
  }

  result += norm(singletAmp)+ norm(triplet0Amp)+ norm(tripletp1Amp)+ norm(tripletm1Amp); 

  }


  if(_usePhasespace){
    result = result + theParamVal.otherParams[_phasespaceKey];
  }
  return result;  

}



void pbarpBaseLh::getDefaultParams(fitParams& fitVal, fitParams& fitErr){ 

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;

  double magFactor=1./sqrt(_jpclsStates.size());  
  std::vector< boost::shared_ptr<const JPCLS> >::iterator it;
  for ( it = _jpclsStates.begin(); it!=_jpclsStates.end(); ++it){
    currentMagValMap[*it] = magFactor;
    currentPhiValMap[*it] = 0.;
    currentMagErrMap[*it] = magFactor;
    currentPhiErrMap[*it] = 0.3;
  }
  fitVal.Mags["pbarp"]=currentMagValMap;
  fitVal.Phis["pbarp"]=currentPhiValMap;
  fitErr.Mags["pbarp"]=currentMagErrMap;
  fitErr.Phis["pbarp"]=currentPhiErrMap;

  std::vector< boost::shared_ptr<AbsXdecAmp> >::iterator itDecs;
  for(itDecs=_decAmps.begin(); itDecs!=_decAmps.end(); ++itDecs){
    (*itDecs)->getDefaultParams(fitVal, fitErr);
  }
}

void pbarpBaseLh::print(std::ostream& os) const{

}


void  pbarpBaseLh::initialize(){

  std::vector<Particle*> fsParticles=pbarpEnv::instance()->finalStateParticles();
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
  
  _pbarpReactionPtr = pbarpEnv::instance()->reaction();
  _jpclsStates=_pbarpReactionPtr->jpclsStates();

}

void pbarpBaseLh::checkParamVariation(fitParams& theParamVal){
  std::vector< boost::shared_ptr<AbsXdecAmp> >::iterator it;
  for (it=_decAmps.begin(); it!=_decAmps.end(); ++it){
    (*it)->checkRecalculation(theParamVal);
  }
}

void pbarpBaseLh::cacheTheAmps(){
  std::vector< boost::shared_ptr<AbsXdecAmp> >::iterator it;
  for (it=_decAmps.begin(); it!=_decAmps.end(); ++it){
    (*it)->cacheAmplitudes();
  }
}

void pbarpBaseLh::updateFitParams(fitParams& theParamVal){

   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > magMap=theParamVal.Mags["pbarp"];
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > phiMap=theParamVal.Phis["pbarp"];

   std::vector< boost::shared_ptr<const JPCLS> >::iterator it;
   for (it=_jpclsStates.begin(); it!=_jpclsStates.end(); ++it){
     double theMag=magMap[*it];
     double thePhi=phiMap[*it];
     _currentParamMags[*it]=theMag;
     _currentParamPhis[*it]=thePhi;
   }

  std::vector< boost::shared_ptr<AbsXdecAmp> >::iterator itDecs;
  for(itDecs=_decAmps.begin(); itDecs!=_decAmps.end(); ++itDecs){
    (*itDecs)->updateFitParams(theParamVal);
  }
}

void pbarpBaseLh::fillMap(std::vector< boost::shared_ptr<const JPCLS> >& pbarpLSs, std::vector<boost::shared_ptr<AbsXdecAmp> >& decs, std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<AbsXdecAmp> >, pawian::Collection::SharedPtrLess >& toFill){


  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  for (itJPCLS = pbarpLSs.begin(); itJPCLS != pbarpLSs.end(); ++itJPCLS){
    std::vector<boost::shared_ptr<AbsXdecAmp> > currentAmpVec;

    std::vector<boost::shared_ptr<AbsXdecAmp> >::iterator itAmp;
    for(itAmp=decs.begin(); itAmp!=decs.end(); ++itAmp){
      if( (*(*itAmp)->jpcPtr())==(*(*itJPCLS)) )  currentAmpVec.push_back(*itAmp);
    }
    toFill[(*itJPCLS)]=currentAmpVec;
  }

}
