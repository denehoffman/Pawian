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

// resBaseLh class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "resUtils/resBaseLh.hh"
#include "resUtils/resReaction.hh"
#include "resUtils/ResChannelEnv.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/FsParticleProjections.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


resBaseLh::resBaseLh(ChannelID channelID) :
  AbsLh(channelID)
  ,_highestJFsp(0)
  ,_isHighestJaPhoton(true)
  ,_withPolarization(false)
{
  initialize();
}


resBaseLh::~resBaseLh()
{;
}


complex<double> resBaseLh::calcSpinDensity(Spin M1, Spin M2, std::string& nameDec, EvtData* theData){

  complex<double> result(0.,0.);

  return result;
}




complex<double> resBaseLh::calcProdPartAmp(Spin lamX, Spin lamDec, std::string nameDec, EvtData* theData,
					     std::map <std::shared_ptr<const JPCLS>,
					     std::vector< std::shared_ptr<AbsXdecAmp> >,
					     pawian::Collection::SharedPtrLess > resAmps){
   complex<double> resultAmp(0.,0.);

   return resultAmp;
}


double resBaseLh::calcEvtIntensity(EvtData* theData, std::shared_ptr<AbsPawianParameters> fitPar){

  double result=0.;

  std::vector< std::shared_ptr<AbsXdecAmp> >::iterator itDecAll;
  for (itDecAll=_decAmps.begin(); itDecAll!=_decAmps.end(); ++itDecAll){
    (*itDecAll)->calcDynamics(theData);
  }


  std::vector< std::vector<Spin> > spinProjections=_fsParticleProjections->spinProjections();
  
  for (unsigned int projId=0; projId<spinProjections.size(); ++projId){
    for (itDecAll=_decAmps.begin(); itDecAll!=_decAmps.end(); ++itDecAll){
      //(*itDecAll)->setSpinProjections(projId);
      std::vector<Spin> currentSpinProjection=spinProjections.at(projId);
      (*itDecAll)->setSpinProjections(currentSpinProjection, projId);
    }

    for( Spin itLam=-_Jmother; itLam<=_Jmother; ++itLam){
      complex<double> lamItAmp(0.,0.);
      
      std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
      for( itDec=_decAmps.begin(); itDec!=_decAmps.end(); ++itDec){
	complex<double> currentDecAmp = (*itDec)->XdecAmp(itLam, theData);
	lamItAmp += currentDecAmp;
      }

      if(_useCohPhasespace) lamItAmp += std::polar( fitPar->Value(_CohPhasespaceKey+"Mag"), fitPar->Value(_CohPhasespaceKey+"Phi") );

      result += fabs(_currentPolVec.at(itLam))*norm(lamItAmp);
    }
  }
  
  if(_usePhasespace) result+=fitPar->Value(_phasespaceKey);

  result *= fitPar->Value(_channelScaleParam);
  return result;

}

void resBaseLh::print(std::ostream& os) const{

}


void  resBaseLh::initialize(){
  const std::shared_ptr<ResChannelEnv> resEnv=std::static_pointer_cast<ResChannelEnv>(GlobalEnv::instance()->ResChannel(_channelID));
  _Jmother = Spin(resEnv->motherParticle()->twoJ(), 2);

  std::vector<Particle*> fsParticles=GlobalEnv::instance()->Channel(_channelID)->finalStateParticles();
  std::vector<Particle*>::iterator itParticle;
  bool highestJFound=false;

  for (itParticle=fsParticles.begin(); itParticle != fsParticles.end(); ++itParticle){
    int current2J = (*itParticle)->twoJ();
    if(current2J>0){
      if(highestJFound){
	Alert << "final states with more than 1 particles with J>0 not supported!!!!" << endmsg;
	exit(1);
      }
      _highestJFsp=Spin(current2J, 2);
      if( (*itParticle)->name() != "photon" ) _isHighestJaPhoton=false;
    }
  }

  _resReactionPtr = std::static_pointer_cast<ResChannelEnv>(GlobalEnv::instance()->ResChannel(_channelID))->reaction();

  // std::vector< std::shared_ptr<IsobarHeliDecay> > theDecs = _resReactionPtr->productionHeliDecays();
  // std::vector< std::shared_ptr<IsobarHeliDecay> >::iterator it;
  // for (it=theDecs.begin(); it!=theDecs.end(); ++it){
  //   //    std::shared_ptr<AbsDecay> currentDec((*it).get() );
  //   //    std::shared_ptr<AbsXdecAmp> currentAmp=XdecAmpRegistry::instance()->getXdecAmp(currentDec);
  //   std::shared_ptr<AbsXdecAmp> currentAmp=XdecAmpRegistry::instance()->getXdecAmp(_channelID, (*it)->absDecPtr());
  //   _decAmps.push_back(currentAmp);
  // }

  std::vector< std::shared_ptr<AbsDecay> > theDecs = resEnv->prodDecayList()->getList();
  std::vector< std::shared_ptr<AbsDecay> >::iterator it;
  for (it=theDecs.begin(); it!=theDecs.end(); ++it){
    std::shared_ptr<AbsXdecAmp> currentAmp=XdecAmpRegistry::instance()->getXdecAmp(_channelID, (*it)->absDecPtr());
    _decAmps.push_back(currentAmp);
  }
  
  _withPolarization=resEnv->polarizedMother();

  int noOfHelStates(2.*_Jmother+1);
  for( Spin itLam=-_Jmother; itLam<=_Jmother; ++itLam){
    std::string currentName=resEnv->motherParticle()->name();
    std::stringstream tmpStrStreamlamJ;
    tmpStrStreamlamJ << itLam;
    currentName+= tmpStrStreamlamJ.str();   
    _polParamNames[itLam]=currentName;
    if(_withPolarization) _currentPolVec[itLam]=1./noOfHelStates;
    else _currentPolVec[itLam]=1.;
  }
}


void resBaseLh::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){

  AbsLh::fillDefaultParams(fitPar);

  if (!_withPolarization) return;

  std::map<Spin, std::string>::iterator it;
  for(it=_polParamNames.begin(); it!=_polParamNames.end(); ++it){
    std::string currentName=it->second;
    double currentVal=_currentPolVec.at(it->first);
    fitPar->Add(currentName, currentVal, 0.05);    
  } 

}

void resBaseLh::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  AbsLh::updateFitParams(fitPar);
  if (!_withPolarization) return;

  std::map<Spin, std::string>::iterator it;
  for(it=_polParamNames.begin(); it!=_polParamNames.end(); ++it){
    std::string currentName=it->second;
    double newVal=fitPar->Value(currentName);
    _currentPolVec.at(it->first)=newVal;    
  }
}


