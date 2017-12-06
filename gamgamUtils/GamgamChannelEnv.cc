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

// GamgamChannelEnv class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "gamgamUtils/GamgamChannelEnv.hh"
#include "ConfigParser/gamgamParser.hh"
#include "gamgamUtils/gamgamReaction.hh"
//#include "epemUtils/epemHist.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/FormationDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/IsobarTensorDecay.hh"
#include "PwaUtils/ProdChannelInfo.hh"
#include "ErrLogger/ErrLogger.hh"





GamgamChannelEnv::GamgamChannelEnv(gamgamParser* theParser) : 
  AbsChannelEnv(theParser, AbsChannelEnv::CHANNEL_GAMGAM)
  ,_theGamGamParser(theParser)
  ,_cmsMass(1.)
  ,_jmax(4)
{
}

void GamgamChannelEnv::setupChannel(ChannelID id){

  AbsChannelEnv::setupGlobal(id);
  
  
  _cmsMass=_theGamGamParser->cmsMass();
  _cmEnergy = _cmsMass;
  _jmax=_theGamGamParser->jMax();
  std::vector<std::string>::const_iterator itStr;
  
  
  //gamgam reaction
  _gamgamReaction=std::shared_ptr<gamgamReaction>(new gamgamReaction(_prodChannelInfoList, id, _jmax));
  std::string dynTypeDefault="WoDynamics";
  
  if (_theGamGamParser->productionFormalism()=="Formation"){
    std::vector< std::shared_ptr<FormationDecay> > prodDecs = _gamgamReaction->formationDecays();
    std::vector< std::shared_ptr<FormationDecay> >::iterator itDec;

    for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
      //      if((*itDec)->prodChannelInfo()->withProdBarrier()) (*itDec)->enableProdBarrier();
      //      else (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
      _prodDecList->addDecay(*itDec);
      InfoMsg <<"added production dacay " << (*itDec)->name() << endmsg;
    }
  }
  else{
    Alert <<"production formalism with the name " << _theGamGamParser->productionFormalism() << " doesn't exist for gam gam reactions!!!" << endmsg;
    exit(1);
  }
  
  //set prefactor for production and decay amplitudes
  AbsChannelEnv::setPrefactors();
  
  //replace suffixes for fit parameter
  AbsChannelEnv::replaceParameterSuffixes();

  //replace mass key
  AbsChannelEnv::replaceMassKeys();

  //add dynamics
  AbsChannelEnv::addDynamics();

  //set decay levels
  AbsChannelEnv::setDecayLevels();
}


