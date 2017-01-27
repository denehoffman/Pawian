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

// ResChannelEnv class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#include "Particle/ParticleTable.hh"
#include "resUtils/ResChannelEnv.hh"
#include "ConfigParser/resParser.hh"
#include "resUtils/resReaction.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/ProdChannelInfo.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"

ResChannelEnv::ResChannelEnv(resParser* theParser) 
  : AbsChannelEnv(theParser, AbsChannelEnv::CHANNEL_RES)
  ,_theResParser(theParser)
  ,_polarizedMother(false)
{
}



void ResChannelEnv::setupChannel(ChannelID id){

   AbsChannelEnv::setupGlobal(id);

   //  // has to be set via parser !!!!
  // double totalyMom=0.04;
  // _initial4Vec = Vector4<double>( sqrt(_cmsMass*_cmsMass+totalyMom*totalyMom), 0., totalyMom, 0.);

  _motherParticle = GlobalEnv::instance()->particleTable()->particle(_theResParser->motherResName());
  if(0==_motherParticle){
    Alert << "mother particle with name\t" << _theResParser->motherResName() << "\tdoesn't exist" << endmsg;
    exit(0);
  }

  _cmEnergy=_motherParticle->mass();

  _polarizedMother=_theResParser->polarizedMother();
  if(_polarizedMother) {
    InfoMsg << "polarization of the mother particle will be taken into accout!!!" << endmsg;
  }

  //  double theMotherMass=_motherParticle->mass();
  _initial4Vec = Vector4<double>(_motherParticle->mass(), 0. ,0. ,0.);
  std::vector<std::string>::const_iterator itStr;


  //epem reaction
  _resReaction=std::shared_ptr<resReaction>(new resReaction(_motherParticle, _prodChannelInfoList, id));

  //fill prodDecayList
  std::vector<std::string> additionalStringVecDummy;
  std::string dynTypeDefault="WoDynamics";

  if(_theResParser->productionFormalism()=="Heli"){
    std::vector< std::shared_ptr<IsobarHeliDecay> > prodDecs= _resReaction->productionHeliDecays();
    std::vector< std::shared_ptr<IsobarHeliDecay> >::iterator itDec;
    for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
      (*itDec)->disableIsospin();
      if((*itDec)->prodChannelInfo()->withProdBarrier()) (*itDec)->enableProdBarrier();
      else (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
      _prodDecList->addDecay(*itDec);
    }
  }
  else if(_theResParser->productionFormalism()=="Cano"){
    std::vector< std::shared_ptr<IsobarLSDecay> > prodDecs= _resReaction->productionCanoDecays();
    std::vector< std::shared_ptr<IsobarLSDecay> >::iterator itDec;
    for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
      (*itDec)->disableIsospin();
      if((*itDec)->prodChannelInfo()->withProdBarrier()) (*itDec)->enableProdBarrier();
      else (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
      _prodDecList->addDecay(*itDec);
    }
  }
  else{
    Alert << "production formalism: " << _theResParser->productionFormalism() << " is not supported for res reactions" << endmsg;
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



