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

// pbarpReaction class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>

#include "pbarpUtils/pbarpReaction.hh"
#include "pbarpUtils/pbarpEnv.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/IsobarTensorDecay.hh"
#include "PwaUtils/ParserBase.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "pbarpUtils/pbarpStatesLS.hh"

pbarpReaction::pbarpReaction(std::vector<std::pair<Particle*, Particle*> >& prodPairs, int lmax) :
  _lmax(lmax)
{
  std::shared_ptr<pbarpStatesLS> thepbarpStates(new pbarpStatesLS(lmax));
  std::vector< std::shared_ptr<const jpcRes> > pbarpJPCStatesAll= thepbarpStates->jpcStates(); 
  std::vector< std::shared_ptr<const jpcRes> >::const_iterator itJPC;
  std::vector< std::shared_ptr<const JPCLS> > all_JPCLSs= thepbarpStates->all_JPCLS_States();
  std::vector< std::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  for(itJPCLS = all_JPCLSs.begin(); itJPCLS != all_JPCLSs.end(); ++itJPCLS){
     std::shared_ptr<const jpcRes> currentJPC = (*itJPCLS);

     if(std::find(_jpcToJPCLSMap[currentJPC].begin(), _jpcToJPCLSMap[currentJPC].end(), (*itJPCLS)) == _jpcToJPCLSMap[currentJPC].end()){
	_jpcToJPCLSMap[currentJPC].push_back(*itJPCLS);
     }

  }

  for(itJPC = pbarpJPCStatesAll.begin(); itJPC!=pbarpJPCStatesAll.end(); ++itJPC){
    bool acceptJPC=false;
    std::vector<std::pair<Particle*, Particle*> >::iterator itPartPairs;
    for (itPartPairs=prodPairs.begin(); itPartPairs!= prodPairs.end(); ++itPartPairs){
      std::string decName=(*itJPC)->name();
      std::shared_ptr<IsobarLSDecay> currentDec(new IsobarLSDecay( (*itJPC),itPartPairs->first, itPartPairs->second, pbarpEnv::instance(), decName));

      if(!currentDec->JPCLSAmps().size()>0)
	 continue;

      bool acceptProd=false;
      for(auto lIt = _jpcToJPCLSMap[*itJPC].begin(); lIt != _jpcToJPCLSMap[*itJPC].end(); ++lIt){
	 if(CheckJPCLSForParticle((std::string&)(*itPartPairs).first->name(), *lIt) &&
	    CheckJPCLSForParticle((std::string&)(*itPartPairs).second->name(), *lIt)){
	    if(std::find(_pbarpJPCLSs.begin(), _pbarpJPCLSs.end(), *lIt) == _pbarpJPCLSs.end())
	       _pbarpJPCLSs.push_back(*lIt);
	    acceptProd = true;
	 }
      }

      if(!acceptProd)
	 continue;

      acceptJPC=true;

      if(pbarpEnv::instance()->parser()->productionFormalism() == "Cano"){
	_prodDecs.push_back(currentDec);
      }
      
      else if(pbarpEnv::instance()->parser()->productionFormalism() == "Tensor"){
	 std::shared_ptr<IsobarTensorDecay> currentTensorDec(new IsobarTensorDecay( (*itJPC),itPartPairs->first, itPartPairs->second, pbarpEnv::instance(), decName));
	 _prodTensorDecs.push_back(currentTensorDec);
      }
      else if(pbarpEnv::instance()->parser()->productionFormalism() == "Heli"){
	 std::shared_ptr<IsobarHeliDecay> currentHeliDec(new IsobarHeliDecay( (*itJPC),itPartPairs->first, itPartPairs->second, pbarpEnv::instance(), decName));
	 _prodHeliDecs.push_back(currentHeliDec);
      }
    }
    if(acceptJPC)
       _pbarpJPCs.push_back(*itJPC);
  }

  std::vector< std::shared_ptr<const JPCLS> > all_pbarpSingletLS = thepbarpStates->singlet_JPCLS_States();
  _pbarpJPCLSsinglet =  extractStates(_pbarpJPCLSs, all_pbarpSingletLS);

  std::vector< std::shared_ptr<const JPCLS> > all_pbarpTriplet0LS = thepbarpStates->triplet0_JPCLS_States();
  _pbarpJPCLStriplet0 =  extractStates(_pbarpJPCLSs, all_pbarpTriplet0LS);

  std::vector< std::shared_ptr<const JPCLS> > all_pbarpTripletp1LS = thepbarpStates->tripletp1_JPCLS_States();
  _pbarpJPCLStripletp1 =  extractStates(_pbarpJPCLSs, all_pbarpTripletp1LS);

  std::vector< std::shared_ptr<const JPCLS> > all_pbarpTripletm1LS = thepbarpStates->tripletm1_JPCLS_States();
  _pbarpJPCLStripletm1 =  extractStates(_pbarpJPCLSs, all_pbarpTripletm1LS);
}



bool pbarpReaction::CheckJPCLSForParticle(std::string& particleName, std::shared_ptr<const JPCLS> theJPCLS){

   std::map<std::string, std::vector<short> > ldropMap = pbarpEnv::instance()->dropPbarpLForParticleData();

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



pbarpReaction::~pbarpReaction(){
}

void pbarpReaction::print(std::ostream& os) const{
  os << "\n pbarp reaction\n";
 
  os << "\n **** JPC states for pbarp system ****** \n";
  std::vector<std::shared_ptr<const jpcRes> >::const_iterator itJPC;
  for(itJPC=_pbarpJPCs.begin(); itJPC!=_pbarpJPCs.end(); ++itJPC){
    (*itJPC)->print(os);
  }

  os << "\n ***** decay chains *******\n";
  if (pbarpEnv::instance()->parser()->productionFormalism() == "Cano"){
    std::vector< std::shared_ptr<IsobarLSDecay> >::const_iterator itIso;
    for( itIso=_prodDecs.begin(); itIso!=_prodDecs.end(); ++itIso){
      (*itIso)->print(os);
    } 
  } 
  else if (pbarpEnv::instance()->parser()->productionFormalism() == "Heli"){
    std::vector< std::shared_ptr<IsobarHeliDecay> >::const_iterator itIso;
    for( itIso=_prodHeliDecs.begin(); itIso!=_prodHeliDecs.end(); ++itIso){
      (*itIso)->print(os);
    } 
  }
  else if (pbarpEnv::instance()->parser()->productionFormalism() == "Tensor"){
    std::vector< std::shared_ptr<IsobarTensorDecay> >::const_iterator itIso;
    for( itIso=_prodTensorDecs.begin(); itIso!=_prodTensorDecs.end(); ++itIso){
      (*itIso)->print(os);
    } 
  }  
}
