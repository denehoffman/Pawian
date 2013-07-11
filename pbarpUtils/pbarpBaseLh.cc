//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
//									  //
//  This file is part of Pawian.					  //
//									  //
//  Pawian is free software: you can redistribute it and/or modify	  //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	  //
//									  //
//  Pawian is distributed in the hope that it will be useful,		  //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  //
//  GNU General Public License for more details.	      		  //
//									  //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	  //
//									  //
//************************************************************************//

// pbarpBaseLh class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "pbarpUtils/pbarpBaseLh.hh"
#include "pbarpUtils/pbarpEnv.hh"
#include "pbarpUtils/pbarpReaction.hh"
#include "pbarpUtils/pbarpEnv.hh"
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


pbarpBaseLh::pbarpBaseLh() :
  AbsLh(pbarpEnv::instance())
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
					     std::map <std::shared_ptr<const JPCLS>,
					     std::vector< std::shared_ptr<AbsXdecAmp> >,
					     pawian::Collection::SharedPtrLess > pbarpAmps){
   complex<double> resultAmp(0.,0.);

   std::map <std::shared_ptr<const JPCLS>,
             std::vector< std::shared_ptr<AbsXdecAmp> >,
             pawian::Collection::SharedPtrLess >::iterator it;

   for(it=pbarpAmps.begin(); it!=pbarpAmps.end(); ++it){
     complex<double> tmpAmp(0.,0.);
     std::shared_ptr<const JPCLS> theJPCLS=it->first;
     std::vector<std::shared_ptr<AbsXdecAmp> > decAmps=it->second;
     std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;

     for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
       Particle* particle1 = (*itDec)->absDec()->daughter1Part();
       Particle* particle2 = (*itDec)->absDec()->daughter2Part();

       if(particle1->name() == nameDec)
	 tmpAmp+=(*itDec)->XdecPartAmp(lamX, lamDec, 1,  theData, 0);
       else if(particle2->name() == nameDec)
	 tmpAmp+=(*itDec)->XdecPartAmp(lamX, lamDec, 2,  theData, 0);
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

  std::map <std::shared_ptr<const JPCLS>, std::vector< std::shared_ptr<AbsXdecAmp> >, pawian::Collection::SharedPtrLess >::iterator it;

  int lamSteps=1;
  if(_isHighestJaPhoton) lamSteps=2;

  for (int lamHigestJFsp=-_highestJFsp; lamHigestJFsp<=_highestJFsp; lamHigestJFsp=lamHigestJFsp+lamSteps){
 
  complex<double> singletAmp(0.,0.);
  for(it=_decAmpsSinglet.begin(); it!=_decAmpsSinglet.end(); ++it){
    complex<double> tmpAmp(0.,0.);
    std::shared_ptr<const JPCLS> theJPCLS=it->first;
    std::vector<std::shared_ptr<AbsXdecAmp> > decAmps=it->second;    
     std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
     for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
       complex<double> currentDecAmp=(*itDec)->XdecAmp(0, theData, lamHigestJFsp);
       std::string currentName=(*itDec)->name();

       std::shared_ptr<AbsDecay> theDec=(*itDec)->absDec();
       // double cgIso0=Clebsch(theDec->iDaughter1(), theDec->i3Daughter1(), theDec->iDaughter2(), theDec->i3Daughter2(), 0, 0);
       // double cgIso1=Clebsch(theDec->iDaughter1(), theDec->i3Daughter1(), theDec->iDaughter2(), theDec->i3Daughter2(), 1, 0);
       // double theIso0=_currentParamIsos[theJPCLS];
       // double theIso1=sqrt(1.-theIso0*theIso0);
       // tmpAmp+=(theIso0*_iso0Map[currentName]+theIso1*_iso1Map[currentName])*currentDecAmp;
       tmpAmp+=_isoMap[theDec->name()]*currentDecAmp;
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
    std::shared_ptr<const JPCLS> theJPCLS=it->first;
    std::vector<std::shared_ptr<AbsXdecAmp> > decAmps=it->second;    
     std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
     for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
       complex<double> currentDecAmp=(*itDec)->XdecAmp(0, theData, lamHigestJFsp);
       std::string currentName=(*itDec)->name();
       
       std::shared_ptr<AbsDecay> theDec=(*itDec)->absDec();
       // double cgIso0=Clebsch(theDec->iDaughter1(), theDec->i3Daughter1(), theDec->iDaughter2(), theDec->i3Daughter2(), 0, 0);
       // double cgIso1=Clebsch(theDec->iDaughter1(), theDec->i3Daughter1(), theDec->iDaughter2(), theDec->i3Daughter2(), 1, 0);
       // double theIso0=_currentParamIsos[theJPCLS];
       // double theIso1=sqrt(1.-theIso0*theIso0);
       //       tmpAmp+=(theIso0*_iso0Map[currentName]+theIso1*_iso1Map[currentName])*currentDecAmp;
       tmpAmp+=_isoMap[theDec->name()]*currentDecAmp;
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
    std::shared_ptr<const JPCLS> theJPCLS=it->first;
    std::vector<std::shared_ptr<AbsXdecAmp> > decAmps=it->second;    
     std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
     for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
       complex<double> currentDecAmp=(*itDec)->XdecAmp(1, theData, lamHigestJFsp);
       std::string currentName=(*itDec)->name();
       std::shared_ptr<AbsDecay> theDec=(*itDec)->absDec();
       // double cgIso0=Clebsch(theDec->iDaughter1(), theDec->i3Daughter1(), theDec->iDaughter2(), theDec->i3Daughter2(), 0, 0);
       // double cgIso1=Clebsch(theDec->iDaughter1(), theDec->i3Daughter1(), theDec->iDaughter2(), theDec->i3Daughter2(), 1, 0);
       // double theIso0=_currentParamIsos[theJPCLS];
       // double theIso1=sqrt(1.-theIso0*theIso0);
       //       tmpAmp+=(theIso0*_iso0Map[currentName]+theIso1*_iso1Map[currentName])*currentDecAmp;
       tmpAmp+=_isoMap[theDec->name()]*currentDecAmp;
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
    std::shared_ptr<const JPCLS> theJPCLS=it->first;
    std::vector<std::shared_ptr<AbsXdecAmp> > decAmps=it->second;    
     std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
     for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
       complex<double> currentDecAmp=(*itDec)->XdecAmp(-1, theData, lamHigestJFsp);
       std::string currentName=(*itDec)->name();
       std::shared_ptr<AbsDecay> theDec=(*itDec)->absDec();
       // double cgIso0=Clebsch(theDec->iDaughter1(), theDec->i3Daughter1(), theDec->iDaughter2(), theDec->i3Daughter2(), 0, 0);
       // double cgIso1=Clebsch(theDec->iDaughter1(), theDec->i3Daughter1(), theDec->iDaughter2(), theDec->i3Daughter2(), 1, 0);
       // double theIso0=_currentParamIsos[theJPCLS];
       // double theIso1=sqrt(1.-theIso0*theIso0);
       //       tmpAmp+=(theIso0*_iso0Map[currentName]+theIso1*_iso1Map[currentName])*currentDecAmp;
       tmpAmp+=_isoMap[theDec->name()]*currentDecAmp;
     }

     double theMag=_currentParamMags[theJPCLS];
     double thePhi=_currentParamPhis[theJPCLS];
     complex<double> expi(cos(thePhi), sin(thePhi));
     tmpAmp*=theJPCLS->preFactor*theMag*expi;
     tripletm1Amp+=tmpAmp;
  }

  result += norm(singletAmp)+ norm(triplet0Amp)+ norm(tripletp1Amp)+ norm(tripletm1Amp); 
  }

  if(_usePhasespace) result+=theParamVal.otherParams[_phasespaceKey];
  return result;  

}



void pbarpBaseLh::getDefaultParams(fitParams& fitVal, fitParams& fitErr){ 

  AbsLh::getDefaultParams(fitVal, fitErr);

  std::map<std::string, std::pair<std::string, std::string> >::iterator itStringMapPair;
  for(itStringMapPair=_iso0Iso1NameMap.begin(); itStringMapPair!=_iso0Iso1NameMap.end(); ++itStringMapPair){
    _isoMap[itStringMapPair->second.first]=sqrt(0.5);
    _isoMap[itStringMapPair->second.second]=sqrt(1.-0.5*0.5);
    fitVal.otherParams["Iso"+itStringMapPair->first+"Range01"]=sqrt(0.5);
    fitErr.otherParams["Iso"+itStringMapPair->first+"Range01"]=0.2;
  }

  std::map< std::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
  std::map< std::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
  std::map< std::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
  std::map< std::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;

  double magFactor=1./sqrt(_jpclsStates.size());  
  std::vector< std::shared_ptr<const JPCLS> >::iterator it;
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
  _igjpcStates=_pbarpReactionPtr->igjpcStates();

}

void pbarpBaseLh::updateFitParams(fitParams& theParamVal){

  AbsLh::updateFitParams(theParamVal);

  std::map<std::string, std::pair<std::string, std::string> >::iterator itStringMapPair;
  for(itStringMapPair=_iso0Iso1NameMap.begin(); itStringMapPair!=_iso0Iso1NameMap.end(); ++itStringMapPair){
    double theVal=theParamVal.otherParams["Iso"+itStringMapPair->first+"Range01"];
    if(theVal<0. || theVal>1.){
      Alert << "fit parameter must be between 0. and 1." << endmsg;
      exit(0);
    }
    _isoMap[itStringMapPair->second.first]=theVal;
    _isoMap[itStringMapPair->second.second]=sqrt(1.-theVal*theVal);
  }
  //  std::map< std::shared_ptr<AbsXdecAmp>, std::shared_ptr<AbsXdecAmp>, pawian::Collection::SharedPtrLess > _iso0Iso1DecAmpMap; 
  // std::map< std::shared_ptr<const jpcRes>, double, pawian::Collection::SharedPtrLess > isoMap=theParamVal.Isos["pbarp"];
  // std::vector< std::shared_ptr<const IGJPC> >::iterator itIGJPC;
  // for (itIGJPC=_igjpcStates.begin(); itIGJPC!=_igjpcStates.end(); ++itIGJPC){
  //   double theVal=isoMap[*itIGJPC];
  //   _currentParamIsos[*itIGJPC]=theVal;
  //  }

   std::map< std::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > magMap=theParamVal.Mags["pbarp"];
   std::map< std::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > phiMap=theParamVal.Phis["pbarp"];
   std::vector< std::shared_ptr<const JPCLS> >::iterator it;
   for (it=_jpclsStates.begin(); it!=_jpclsStates.end(); ++it){
     double theMag=magMap[*it];
     double thePhi=phiMap[*it];
     _currentParamMags[*it]=theMag;
     _currentParamPhis[*it]=thePhi;
   }

}

void pbarpBaseLh::fillMap(std::vector< std::shared_ptr<const JPCLS> >& pbarpLSs, std::vector<std::shared_ptr<AbsXdecAmp> >& decs, std::map< std::shared_ptr<const JPCLS>, std::vector<std::shared_ptr<AbsXdecAmp> >, pawian::Collection::SharedPtrLess >& toFill){


  std::vector< std::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  for (itJPCLS = pbarpLSs.begin(); itJPCLS != pbarpLSs.end(); ++itJPCLS){
    std::vector<std::shared_ptr<AbsXdecAmp> > currentAmpVec;

    std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itAmp;
    for(itAmp=decs.begin(); itAmp!=decs.end(); ++itAmp){
      if( (*(*itAmp)->jpcPtr())==(*(*itJPCLS)) )  currentAmpVec.push_back(*itAmp);
    }
    toFill[(*itJPCLS)]=currentAmpVec;
  }

}

void pbarpBaseLh::fillIsos(){
  //first look for iso0 decay amplitudes and fir it in the map
  std::vector< std::shared_ptr<AbsXdecAmp> >::iterator it;
  for(it= _decAmps.begin(); it!= _decAmps.end(); ++it){
    if((*it)->absDec()->motherIGJPC()->I==0){
      _iso0DecAmps.push_back(*it);
      _isoMap[(*it)->name()]=1.;
    }
    if((*it)->absDec()->motherIGJPC()->I==1){
      _iso1DecAmps.push_back(*it);
      _isoMap[(*it)->name()]=1.;
    }
  }

  //now look if iso0 and iso1 exisist 
  for(it= _iso0DecAmps.begin(); it!= _iso0DecAmps.end(); ++it){
    std::vector< std::shared_ptr<AbsXdecAmp> >::iterator it2;
    for(it2= _iso1DecAmps.begin(); it2!= _iso1DecAmps.end(); ++it2){
      if((*it)->jpcDecsName()==(*it2)->jpcDecsName()){
	//found
	std::pair<std::string, std::string> currentNamePair=make_pair((*it)->name(), (*it2)->name());
	_iso0Iso1NameMap[(*it)->isoKey()]=currentNamePair;
      }
    }
  } 
}
