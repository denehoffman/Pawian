//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// gamgamBaseLh class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "gamgamUtils/gamgamBaseLh.hh"
#include "gamgamUtils/gamgamReaction.hh"
#include "gamgamUtils/GamgamChannelEnv.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "PwaUtils/FsParticleProjections.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

gamgamBaseLh::gamgamBaseLh(ChannelID channelID) :
  AbsLh(channelID)
  ,_gamgamChannelEnv(std::static_pointer_cast<GamgamChannelEnv> (GlobalEnv::instance()->GamgamChannel(channelID)))
{
  initialize();
}


gamgamBaseLh::~gamgamBaseLh()
{;
}


complex<double> gamgamBaseLh::calcSpinDensity(Spin M1, Spin M2, std::string& nameDec, EvtData* theData){

  complex<double> result(0.,0.);

  return result;
}




complex<double> gamgamBaseLh::calcProdPartAmp(Spin lamX, Spin lamDec, std::string nameDec, EvtData* theData,
					     std::map <std::shared_ptr<const JPCLS>,
					     std::vector< std::shared_ptr<AbsXdecAmp> >,
					     pawian::Collection::SharedPtrLess > pbarpAmps){
   complex<double> resultAmp(0.,0.);

   return resultAmp;
}


double gamgamBaseLh::calcEvtIntensity( EvtData* theData, std::shared_ptr<AbsPawianParameters> fitPar){

  double result=0.;

  //  InfoMsg << "_decAmps.size(): " << _decAmps.size() << endmsg;

  std::vector< std::shared_ptr<AbsXdecAmp> >::iterator itDecAll;
  for (itDecAll=_decAmps.begin(); itDecAll!=_decAmps.end(); ++itDecAll){
    (*itDecAll)->calcDynamics(theData);
  }

  std::vector< std::vector<Spin> > spinProjections=_fsParticleProjections->spinProjections();
  
  for (unsigned int projId=0; projId<spinProjections.size(); ++projId){
    for (itDecAll=_decAmps.begin(); itDecAll!=_decAmps.end(); ++itDecAll){
      (*itDecAll)->setSpinProjections(projId);
    }
    
    complex<double> lamX0Amp(0.,0.);
    
    std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
    Spin lamX=0;
    for( itDec=_decAmps.begin(); itDec!=_decAmps.end(); ++itDec){
      complex<double> currentDecAmp=(*itDec)->XdecAmp(lamX, theData);
      //      InfoMsg << "currentDecAmpX0: " << currentDecAmp << endmsg;
      lamX0Amp+=currentDecAmp;
    }


    complex<double> lamXp2Amp(0.,0.);
    lamX=2;
    for( itDec=_decAmps.begin(); itDec!=_decAmps.end(); ++itDec){
      complex<double> currentDecAmp=(*itDec)->XdecAmp(lamX, theData);
      //      InfoMsg << "currentDecAmpX2: " << currentDecAmp << endmsg;
      lamXp2Amp+=currentDecAmp;
    }

    complex<double> lamXm2Amp(0.,0.);
    lamX=-2;
    for( itDec=_decAmps.begin(); itDec!=_decAmps.end(); ++itDec){
      complex<double> currentDecAmp=(*itDec)->XdecAmp(lamX, theData);
      //      InfoMsg << "currentDecAmpX-2: " << currentDecAmp << endmsg;
      lamXm2Amp+=currentDecAmp;
    }

   
    if(_useCohPhasespace){
      lamX0Amp += std::polar( fitPar->Value(_CohPhasespaceKey+"Mag"), fitPar->Value(_CohPhasespaceKey+"Phi") );
      lamXp2Amp += std::polar( fitPar->Value(_CohPhasespaceKey+"Mag"), fitPar->Value(_CohPhasespaceKey+"Phi") );
      lamXm2Amp += std::polar( fitPar->Value(_CohPhasespaceKey+"Mag"), fitPar->Value(_CohPhasespaceKey+"Phi") );
    }
 
    result += norm(lamX0Amp) + norm(lamXp2Amp) + norm(lamXm2Amp);
  }
  
  if(_usePhasespace) result+=fitPar->Value(_phasespaceKey);

  result *= fitPar->Value(_channelScaleParam);
  //  InfoMsg << "result: " << result << endmsg;
  return result;
}

void gamgamBaseLh::print(std::ostream& os) const{

}


void  gamgamBaseLh::initialize(){
  std::vector< std::shared_ptr<AbsDecay> > theDecs = _gamgamChannelEnv->prodDecayList()->getList();
  std::vector< std::shared_ptr<AbsDecay> >::iterator it;
  for (it=theDecs.begin(); it!=theDecs.end(); ++it){
    InfoMsg << "theDecs->name: " << (*it)->name() << endmsg;
    std::shared_ptr<AbsXdecAmp> currentAmp=XdecAmpRegistry::instance()->getXdecAmp(_channelID, (*it)->absDecPtr());
    _decAmps.push_back(currentAmp);
  }
}




