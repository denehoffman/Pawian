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

// EpemChannelEnv class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "epemUtils/EpemChannelEnv.hh"
#include "ConfigParser/epemParser.hh"
#include "epemUtils/epemReaction.hh"
//#include "epemUtils/epemHist.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/IsobarTensorDecay.hh"
#include "PwaUtils/ProdChannelInfo.hh"
#include "ErrLogger/ErrLogger.hh"





EpemChannelEnv::EpemChannelEnv(epemParser* theParser) : 
  AbsChannelEnv(theParser, AbsChannelEnv::CHANNEL_EPEM)
  ,_theEpEmParser(theParser)
  ,_cmsMass(1.)
{
}

void EpemChannelEnv::setupChannel(ChannelID id){

  AbsChannelEnv::setupGlobal(id);
  
  
  _cmsMass=_theEpEmParser->cmsMass();
  _cmEnergy = _cmsMass;
  // has to be set via parser !!!!
  double totalyMom=0.04;
  _initial4Vec = Vector4<double>( sqrt(_cmsMass*_cmsMass+totalyMom*totalyMom), 0., totalyMom, 0.);
  
  std::vector<std::string>::const_iterator itStr;
  
  
  //epem reaction
  _epemReaction=std::shared_ptr<epemReaction>(new epemReaction(_prodChannelInfoList, id));
  std::vector<std::string> additionalStringVecDummy;
  std::string dynTypeDefault="WoDynamics";
  
  //  std::vector< std::shared_ptr<AbsDecay> > prodDecs;
  if (_theEpEmParser->productionFormalism()=="Heli" || _theEpEmParser->productionFormalism()=="HeliMultipole"){
    std::vector< std::shared_ptr<IsobarHeliDecay> > prodDecs;
    if (_theEpEmParser->productionFormalism()=="Heli") prodDecs = _epemReaction->productionHeliDecays();
    else prodDecs = _epemReaction->productionHeliMultipoleDecays();
    std::vector< std::shared_ptr<IsobarHeliDecay> >::iterator itDec;
    for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
      if((*itDec)->prodChannelInfo()->withProdBarrier()) (*itDec)->enableProdBarrier();
      else (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
      _prodDecList->addDecay(*itDec);
    }
  }
  else if (_theEpEmParser->productionFormalism()=="Tensor"){
    std::vector< std::shared_ptr<IsobarTensorDecay> > prodDecs;
    if(_theEpEmParser->productionTensorRadType() == "Zou"){
      InfoMsg << "use productionTensorRadType " << _theEpEmParser->productionTensorRadType() << endmsg;
      prodDecs=_epemReaction->productionTensorZouDecays();
    } 
    else if(_theEpEmParser->productionTensorRadType() == "Default"){
      InfoMsg << "use productionTensorRadType " << _theEpEmParser->productionTensorRadType() << endmsg;
      prodDecs=_epemReaction->productionTensorDecays(); //default
    } 
    else{
      Alert <<"productionTensorRadType with the name " << _theEpEmParser->productionTensorRadType() << " doesn't exist!!!" << endmsg;
      exit(0);
    }
    
    
    std::vector< std::shared_ptr<IsobarTensorDecay> >::iterator itDec;
    for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
      if((*itDec)->prodChannelInfo()->withProdBarrier()) (*itDec)->enableProdBarrier();
      else (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
      _prodDecList->addDecay(*itDec);
    }
  }
  else if (_theEpEmParser->productionFormalism()=="Cano"){
    std::vector< std::shared_ptr<IsobarLSDecay> > prodDecs;
    if (_theEpEmParser->productionFormalism()=="Cano") prodDecs = _epemReaction->productionCanoDecays();
    
    std::vector< std::shared_ptr<IsobarLSDecay> >::iterator itDec;
    for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
      if((*itDec)->prodChannelInfo()->withProdBarrier()) (*itDec)->enableProdBarrier();
      else (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
      _prodDecList->addDecay(*itDec);
    }
  }
  else{
    Alert <<"production formalism\t" << _theEpEmParser->productionFormalism() << "\t is not supported!!!" << endmsg;
    exit(0);
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


