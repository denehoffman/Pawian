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

// gammapReaction class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>

#include "gammapUtils/gammapReaction.hh"
#include "gammapUtils/GammapChannelEnv.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/IsobarTensorDecay.hh"
#include "ConfigParser/ParserBase.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "gammapUtils/gammapStates.hh"

gammapReaction::gammapReaction(std::vector<std::pair<Particle*, Particle*> >& prodPairs, ChannelID channelID, int lmax) :
   _channelID(channelID)
  ,_lmax(lmax)
{
  std::shared_ptr<gammapStates> thegammapStates(new gammapStates(lmax));
  std::vector< std::shared_ptr<const IGJPC> > gammapIGJPCStatesAll= thegammapStates->igjpcStates();
  std::vector< std::shared_ptr<const IGJPC> >::const_iterator itIGJPC;
  std::vector< std::shared_ptr<const JPCLS> > all_JPCljs= thegammapStates->jpcljStates();
  
  std::vector< std::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::map< std::shared_ptr<const jpcRes>, std::vector< std::shared_ptr<const JPCLS> >, pawian::Collection::SharedPtrLess> currentjpcToJPCljMap;
  for(itJPCLS = all_JPCljs.begin(); itJPCLS != all_JPCljs.end(); ++itJPCLS){
      std::shared_ptr<const jpcRes> currentJPC = (*itJPCLS);

      if(std::find(currentjpcToJPCljMap[currentJPC].begin(), currentjpcToJPCljMap[currentJPC].end(), (*itJPCLS)) == currentjpcToJPCljMap[currentJPC].end()){
   	currentjpcToJPCljMap[currentJPC].push_back(*itJPCLS);
      }

  }

  for(itIGJPC = gammapIGJPCStatesAll.begin(); itIGJPC!=gammapIGJPCStatesAll.end(); ++itIGJPC){
    std::shared_ptr<const jpcRes> currentJPC(new jpcRes((*itIGJPC)));
    bool acceptJPC=false;
    std::vector<std::pair<Particle*, Particle*> >::iterator itPartPairs;
    for (itPartPairs=prodPairs.begin(); itPartPairs!= prodPairs.end(); ++itPartPairs){
      std::string decName=(*itIGJPC)->name();
      std::shared_ptr<IsobarLSDecay> currentDec(new IsobarLSDecay( (*itIGJPC),itPartPairs->first, itPartPairs->second, _channelID, decName));
      currentDec->setProductionAmp();
      currentDec->extractStates();

      if(!currentDec->JPCLSAmps().size()>0 || fabs(currentDec->isospinCG())<1.e-10){
 	InfoMsg << "Skip amplitude " << currentDec->name() << "\tJPCLSAmps().size(): " << currentDec->JPCLSAmps().size() <<"\tisospinCG(): " << currentDec->isospinCG() << endmsg;
  	continue;
      }

      bool acceptProd=false;
      for(auto lIt = currentjpcToJPCljMap[*itIGJPC].begin(); lIt != currentjpcToJPCljMap[*itIGJPC].end(); ++lIt){
  	if(CheckJPCLSForParticle((std::string&)(*itPartPairs).first->name(), *lIt) &&
   	   CheckJPCLSForParticle((std::string&)(*itPartPairs).second->name(), *lIt)){
   	  if(std::find(_gammapJPCljs.begin(), _gammapJPCljs.end(), *lIt) == _gammapJPCljs.end())
   	    _gammapJPCljs.push_back(*lIt);
	  _jpcToJPCljMap[currentJPC].push_back(*lIt);
   	  acceptProd = true;
   	}
      }

      if(!acceptProd)	 continue;

      acceptJPC=true;

      if(GlobalEnv::instance()->parser()->productionFormalism() == "Cano"){
  	_prodLSDecs.push_back(currentDec);
      }
      else if(GlobalEnv::instance()->parser()->productionFormalism() == "Tensor"){
  	std::shared_ptr<IsobarTensorDecay> currentTensorDec(new IsobarTensorDecay( (*itIGJPC),itPartPairs->first, itPartPairs->second, _channelID, decName));
	currentTensorDec->setProductionAmp();
  	currentTensorDec->extractStates();
  	_prodTensorDecs.push_back(currentTensorDec);
      }
      else if(GlobalEnv::instance()->parser()->productionFormalism() == "Heli"){
  	std::shared_ptr<IsobarHeliDecay> currentHeliDec(new IsobarHeliDecay( (*itIGJPC),itPartPairs->first, itPartPairs->second, _channelID, decName));
	currentHeliDec->setProductionAmp();
  	currentHeliDec->extractStates();
  	_prodHeliDecs.push_back(currentHeliDec);
      }
    }
    if(acceptJPC){
      //      std::shared_ptr<const jpcRes> acceptedJPC = (*itIGJPC);
      _gammapIGJPCs.push_back(*itIGJPC);
      _jpcToIGJPCMap[currentJPC].push_back(*itIGJPC);
    }
  }

  // std::vector< std::shared_ptr<const JPCLS> > all_pbarpSingletLS = thepbarpStates->singlet_JPCLS_States();
  // _pbarpJPCLSsinglet =  extractStates(_pbarpJPCLSs, all_pbarpSingletLS);

  // std::vector< std::shared_ptr<const JPCLS> > all_pbarpTriplet0LS = thepbarpStates->triplet0_JPCLS_States();
  // _pbarpJPCLStriplet0 =  extractStates(_pbarpJPCLSs, all_pbarpTriplet0LS);

  // std::vector< std::shared_ptr<const JPCLS> > all_pbarpTripletp1LS = thepbarpStates->tripletp1_JPCLS_States();
  // _pbarpJPCLStripletp1 =  extractStates(_pbarpJPCLSs, all_pbarpTripletp1LS);

  // std::vector< std::shared_ptr<const JPCLS> > all_pbarpTripletm1LS = thepbarpStates->tripletm1_JPCLS_States();
  // _pbarpJPCLStripletm1 =  extractStates(_pbarpJPCLSs, all_pbarpTripletm1LS);
}



bool gammapReaction::CheckJPCLSForParticle(std::string& particleName, std::shared_ptr<const JPCLS> theJPCLS){

   std::map<std::string, std::vector<short> > ldropMap = std::static_pointer_cast<GammapChannelEnv>(GlobalEnv::instance()->Channel(_channelID))->dropGammapLForParticleData();

   auto ldropIt = ldropMap.find(particleName);

   if(ldropIt == ldropMap.end())
       return true;

   for(auto ldropListIt = (*ldropIt).second.begin(); ldropListIt != (*ldropIt).second.end(); ++ldropListIt){
      if((*ldropListIt) == theJPCLS->L){
   	 return false;
      }
   }

   return true;
}



gammapReaction::~gammapReaction(){
}

void gammapReaction::print(std::ostream& os) const{
  os << "\n gamma p reaction\n";

  os << "\n **** IG JPC states for the gamma p system ****** \n";
  std::vector<std::shared_ptr<const IGJPC> >::const_iterator itIGJPC;
  for(itIGJPC=_gammapIGJPCs.begin(); itIGJPC!=_gammapIGJPCs.end(); ++itIGJPC){
    (*itIGJPC)->print(os);
  }

  os << "\n\n**** map JPC to IGJPC ****\n";
  for( std::map< std::shared_ptr<const jpcRes>, std::vector< std::shared_ptr<const IGJPC> >, pawian::Collection::SharedPtrLess>::const_iterator it=_jpcToIGJPCMap.begin(); it!=_jpcToIGJPCMap.end(); ++it){
    os << "\n JPC:";
    it->first->print(os);

    os << "\n IGJPCs:";
    for(std::vector< std::shared_ptr<const IGJPC> >::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2){
      os << "\t";
      (*it2)->print(os);
      os << "\n";
    }
  } 

  os << "\n\n**** map JPC to JPClj ****\n";
  for( std::map< std::shared_ptr<const jpcRes>, std::vector< std::shared_ptr<const JPCLS> >, pawian::Collection::SharedPtrLess>::const_iterator it=_jpcToJPCljMap.begin(); it!=_jpcToJPCljMap.end(); ++it){
    os << "\n JPC:";
    it->first->print(os);

    os << "\n JPClj:";
    for(std::vector< std::shared_ptr<const JPCLS> >::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2){
      os << "\t";
      (*it2)->print(os);
      os << "\n";
    }
  } 

  os << "\n ***** decay chains *******\n";
  if (GlobalEnv::instance()->parser()->productionFormalism() == "Cano"){
    std::vector< std::shared_ptr<IsobarLSDecay> >::const_iterator itIso;
    for( itIso=_prodLSDecs.begin(); itIso!=_prodLSDecs.end(); ++itIso){
      (*itIso)->print(os);
    }
  }
  else if (GlobalEnv::instance()->parser()->productionFormalism() == "Heli"){
    std::vector< std::shared_ptr<IsobarHeliDecay> >::const_iterator itIso;
    for( itIso=_prodHeliDecs.begin(); itIso!=_prodHeliDecs.end(); ++itIso){
      (*itIso)->print(os);
    }
  }
  else if (GlobalEnv::instance()->parser()->productionFormalism() == "Tensor"){
    std::vector< std::shared_ptr<IsobarTensorDecay> >::const_iterator itIso;
    for( itIso=_prodTensorDecs.begin(); itIso!=_prodTensorDecs.end(); ++itIso){
      (*itIso)->print(os);
    }
  }
}
