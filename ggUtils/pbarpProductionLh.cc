//************************************************************************//
//									  //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// pbarProductionLh class definition file. -*- C++ -*-
// Copyright 2020 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "ggUtils/pbarpProductionLh.hh"
#include "ggUtils/ggReaction.hh"
#include "ggUtils/GGChannelEnv.hh"
#include "PwaUtils/pbarpState.hh"
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
#include "ConfigParser/ParserBase.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

pbarpProductionLh::pbarpProductionLh(ChannelID channelID) :
   AbsLh(channelID)
   ,_ggChannelEnv(std::static_pointer_cast<GGChannelEnv> (GlobalEnv::instance()->GGChannel(channelID)))
   ,_useProdDynamics(false) {
   initialize();
}

pbarpProductionLh::~pbarpProductionLh()
{;
}


complex<double> pbarpProductionLh::calcSpinDensity(Spin M1, Spin M2, std::string& nameDec, 
						   EvtData* theData){
  // dummy so far
   complex<double> result(0.,0.);
   return result;
}


complex<double> pbarpProductionLh::calcProdPartAmp(Spin lamX, Spin lamDec, std::string nameDec, 
						   EvtData* theData,
      std::map <std::shared_ptr<const JPCLS>,
      std::vector< std::shared_ptr<AbsXdecAmp> >,
      pawian::Collection::SharedPtrLess > pbarpAmps) {
  // dummy so far
   complex<double> resultAmp(0.,0.);

   return resultAmp;
}


double pbarpProductionLh::calcEvtIntensity( EvtData* theData, std::shared_ptr<AbsPawianParameters> fitPar){

   double result=0.;

   std::vector< std::shared_ptr<AbsXdecAmp> >::iterator itDecAll;
   for (itDecAll=_decAmps.begin(); itDecAll!=_decAmps.end(); ++itDecAll){
      (*itDecAll)->calcDynamics(theData);
   }

   std::vector< std::vector<Spin> > spinProjections=_fsParticleProjections->spinProjections();

   // BEGIN Outer incoherent sum over spinProjections of FS particles...
   for (unsigned int projId=0; projId<spinProjections.size(); ++projId) {
     for (itDecAll=_decAmps.begin(); itDecAll!=_decAmps.end(); ++itDecAll) {
       //(*itDecAll)->setSpinProjections(projId);
       std::vector<Spin> currentSpinProjection=spinProjections.at(projId);
      (*itDecAll)->setSpinProjections(currentSpinProjection, projId);
     }

    complex<double> singletTerm(0.,0.);
    complex<double> triplet0Term(0.,0.);
    complex<double> tripletm1Term(0.,0.);
    complex<double> tripletp1Term(0.,0.);
   Spin lampbarp=0;

    std::vector< std::shared_ptr<AbsXdecAmp> >::iterator it;
    //loop over singlet states
    for(it=_singletAmps.begin(); it!=_singletAmps.end(); ++it){
      singletTerm+=(*it)->XdecAmp(lampbarp, theData);
    }

    //loop over triplet0 state
    for(it=_triplet0Amps.begin(); it!=_triplet0Amps.end(); ++it){
      triplet0Term+=(*it)->XdecAmp(lampbarp, theData);
    }

    //loop over triplet-1 state
    lampbarp=-1;
    for(it=_triplet1Amps.begin(); it!=_triplet1Amps.end(); ++it){
      tripletm1Term+=(*it)->XdecAmp(lampbarp, theData);
    }

    //loop over triplet+1 state
    lampbarp=+1;
    for(it=_triplet1Amps.begin(); it!=_triplet1Amps.end(); ++it){
      tripletp1Term+=(*it)->XdecAmp(lampbarp, theData);
    }

    if(_useCohPhasespace){
      singletTerm += std::polar( fitPar->Value(_CohPhasespaceKey+"Mag"), 
                              fitPar->Value(_CohPhasespaceKey+"Phi") );
      triplet0Term += std::polar( fitPar->Value(_CohPhasespaceKey+"Mag"), 
                               fitPar->Value(_CohPhasespaceKey+"Phi") );
      tripletm1Term += std::polar( fitPar->Value(_CohPhasespaceKey+"Mag"), 
                               fitPar->Value(_CohPhasespaceKey+"Phi") );
      tripletp1Term += std::polar( fitPar->Value(_CohPhasespaceKey+"Mag"), 
                               fitPar->Value(_CohPhasespaceKey+"Phi") );
    }

    result += norm(singletTerm)+norm(triplet0Term)+norm(tripletm1Term)+norm(tripletp1Term);
   }
      
   if(_usePhasespace) result+=fitPar->Value(_phasespaceKey);
   if(_useProdDynamics) result*=fabs(_dyn->eval(theData,0).real());
   result *= fitPar->Value(_channelScaleParam);
   return result;
}

void pbarpProductionLh::print(std::ostream& os) const{

}

void pbarpProductionLh::initialize(){
   std::vector< std::shared_ptr<AbsDecay> > theDecs = _ggChannelEnv->prodDecayList()->getList();
   std::vector< std::shared_ptr<AbsDecay> >::iterator it;
   for (it=theDecs.begin(); it!=theDecs.end(); ++it){
     InfoMsg << "theDecs->name: " << (*it)->name() << endmsg;
     std::shared_ptr<AbsXdecAmp> currentAmp=XdecAmpRegistry::instance()->getXdecAmp(_channelID, 
										    (*it)->absDecPtr());

     //fill pbarp singlet triplet0 and triplet1 lists
     pbarpState thepbarpState(currentAmp->jpcPtr());
     thepbarpState.print(std::cout);

     if( thepbarpState.isSingletState() ) _singletAmps.push_back(currentAmp);
     else if(thepbarpState.isTripletState()){
       if(thepbarpState.isTriplet0State()) _triplet0Amps.push_back(currentAmp);
       if(thepbarpState.isTriplet1State()) _triplet1Amps.push_back(currentAmp);
     }
     else {
       ErrMsg << "pbar p system cannot exhibit the spin quantum numbers: " << currentAmp->jpcPtr()->name() << " " << endmsg;
       exit(1);
     }

     _decAmps.push_back(currentAmp);      
   }

   std::string dynString=_ggChannelEnv->parser()->productionDynamics();
   if(dynString=="woFormPol" || dynString=="FormPol0" || dynString=="FormPol1" || dynString=="FormPol2"){
     _useProdDynamics=true;
     std::shared_ptr<FormationDecay> motherFormDec = _ggChannelEnv->reaction()->motherProdDec();
     std::vector<std::string> additionalStringVecDummy;
     motherFormDec->enableDynamics(dynString, additionalStringVecDummy);
     _dyn=motherFormDec->getDynamics();
   }
}

void pbarpProductionLh::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  if(_useProdDynamics) _dyn->fillDefaultParams(fitPar);
  AbsLh::fillDefaultParams(fitPar);
}

void pbarpProductionLh::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  if(_useProdDynamics) _dyn->updateFitParams(fitPar);
  AbsLh::updateFitParams(fitPar);
}
