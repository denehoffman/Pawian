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

// epemBaseLh class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "epemUtils/epemBaseLh.hh"
#include "epemUtils/epemReaction.hh"
#include "epemUtils/EpemChannelEnv.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

epemBaseLh::epemBaseLh(ChannelID channelID) :
  AbsLh(channelID)
  ,_epemChannelEnv(std::static_pointer_cast<EpemChannelEnv> (GlobalEnv::instance()->EpemChannel(channelID)))
  ,_highestJFsp(0)
  ,_isHighestJaPhoton(true)
{
}


epemBaseLh::~epemBaseLh()
{;
}


complex<double> epemBaseLh::calcSpinDensity(Spin M1, Spin M2, std::string& nameDec, EvtData* theData){

  complex<double> result(0.,0.);

  return result;
}




complex<double> epemBaseLh::calcProdPartAmp(Spin lamX, Spin lamDec, std::string nameDec, EvtData* theData,
					     std::map <std::shared_ptr<const JPCLS>,
					     std::vector< std::shared_ptr<AbsXdecAmp> >,
					     pawian::Collection::SharedPtrLess > pbarpAmps){
   complex<double> resultAmp(0.,0.);

   return resultAmp;
}


// double epemBaseLh::calcEvtIntensity(EvtData* theData, fitParCol& theParamVal){

//   double result=0.;

//   std::vector< std::shared_ptr<AbsXdecAmp> >::iterator itDecAll;
//   for (itDecAll=_decAmps.begin(); itDecAll!=_decAmps.end(); ++itDecAll){
//     (*itDecAll)->calcDynamics(theData);
//   }

//   Spin lamSteps=1;
//   if(_isHighestJaPhoton) lamSteps=2;

//   for (Spin lamHigestJFsp=-_highestJFsp; lamHigestJFsp<=_highestJFsp; lamHigestJFsp=lamHigestJFsp+lamSteps){

//     complex<double> lamp1Amp(0.,0.);

//     std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
//     Spin lamepem=1;
//     for( itDec=_decAmps.begin(); itDec!=_decAmps.end(); ++itDec){
//       complex<double> currentDecAmp=(*itDec)->XdecAmp(lamepem, theData, lamHigestJFsp);
//       lamp1Amp+=currentDecAmp;
//     }


//     complex<double> lamm1Amp(0.,0.);
//     lamepem=-1;
//     for( itDec=_decAmps.begin(); itDec!=_decAmps.end(); ++itDec){
//       complex<double> currentDecAmp=(*itDec)->XdecAmp(lamepem, theData, lamHigestJFsp);
//       lamm1Amp+=currentDecAmp;
//     }

//     result += norm(lamp1Amp) + norm(lamm1Amp);
//   }

//   if(_usePhasespace) result+=theParamVal.otherParams[_phasespaceKey];

//   result *= theParamVal.otherParams.at(_channelScaleParam);

//   return result;

// }

double epemBaseLh::calcEvtIntensity( EvtData* theData, std::shared_ptr<AbsPawianParameters> fitPar){

  double result=0.;

  std::vector< std::shared_ptr<AbsXdecAmp> >::iterator itDecAll;
  for (itDecAll=_decAmps.begin(); itDecAll!=_decAmps.end(); ++itDecAll){
    (*itDecAll)->calcDynamics(theData);
  }

  Spin lamSteps=1;
  if(_isHighestJaPhoton) lamSteps=2;

  for (Spin lamHigestJFsp=-_highestJFsp; lamHigestJFsp<=_highestJFsp; lamHigestJFsp=lamHigestJFsp+lamSteps){

    complex<double> lamp1Amp(0.,0.);

    std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
    Spin lamepem=1;
    for( itDec=_decAmps.begin(); itDec!=_decAmps.end(); ++itDec){
      complex<double> currentDecAmp=(*itDec)->XdecAmp(lamepem, theData, lamHigestJFsp);
      lamp1Amp+=currentDecAmp;
    }


    complex<double> lamm1Amp(0.,0.);
    lamepem=-1;
    for( itDec=_decAmps.begin(); itDec!=_decAmps.end(); ++itDec){
      complex<double> currentDecAmp=(*itDec)->XdecAmp(lamepem, theData, lamHigestJFsp);
      lamm1Amp+=currentDecAmp;
    }

    result += norm(lamp1Amp) + norm(lamm1Amp);
  }

  if(_usePhasespace) result+=fitPar->Value(_phasespaceKey);

  result *= fitPar->Value(_channelScaleParam);

  return result;
}

void epemBaseLh::print(std::ostream& os) const{

}


void  epemBaseLh::initialize(){

  std::vector<Particle*> fsParticles=_epemChannelEnv->finalStateParticles();
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

  //  _epemReactionPtr =  std::static_pointer_cast<EpemChannelEnv>(GlobalEnv::instance()->EpemChannel(_channelID))->reaction();

  std::vector< std::shared_ptr<AbsDecay> > theDecs = _epemChannelEnv->prodDecayList()->getList();
  std::vector< std::shared_ptr<AbsDecay> >::iterator it;
  for (it=theDecs.begin(); it!=theDecs.end(); ++it){
    std::shared_ptr<AbsXdecAmp> currentAmp=XdecAmpRegistry::instance()->getXdecAmp(_channelID, (*it)->absDecPtr());
    _decAmps.push_back(currentAmp);
  }
}




