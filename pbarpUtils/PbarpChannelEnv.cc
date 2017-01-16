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

// PbarpChannelEnv class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "pbarpUtils/PbarpChannelEnv.hh"
#include "ConfigParser/pbarpParser.hh"
#include "pbarpUtils/pbarpReaction.hh"
//#include "pbarpUtils/pbarpHist.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/IsobarTensorDecay.hh"
#include "PwaUtils/ProdChannelInfo.hh"
#include "ErrLogger/ErrLogger.hh"





PbarpChannelEnv::PbarpChannelEnv(pbarpParser* theParser) : AbsChannelEnv(theParser, AbsChannelEnv::CHANNEL_PBARP)
  ,_lmax(0)
  ,_pbarMomentum(0)
  ,_thePbarpParser(theParser)
{
}

void PbarpChannelEnv::setup(ChannelID id){

   AbsChannelEnv::setup(id);

   //Antiproton momentum
  _pbarMomentum = _thePbarpParser->getpbarMomentum();

   double pMass=GlobalEnv::instance()->particleTable()->particle("proton")->mass();
   double antipMass=GlobalEnv::instance()->particleTable()->particle("antiproton")->mass();
   _initial4Vec = Vector4<double>(pMass+sqrt(antipMass*antipMass+_pbarMomentum*_pbarMomentum), 0., 0., _pbarMomentum);
   _cmEnergy = _initial4Vec.M();
     
   //Lmax
   _lmax=_thePbarpParser->getLMax();

   // individual Lmax settings
   std::vector<std::string> theDropPbarpLForParticles = _thePbarpParser->dropPbarpLForParticle();
   for(auto ldropIt = theDropPbarpLForParticles.begin(); ldropIt != theDropPbarpLForParticles.end(); ++ldropIt){
      std::string particle;
      short l;
      std::stringstream stringStr;
      stringStr << (*ldropIt);
      stringStr >> particle >> l;
      _dropPbarpLForParticleData[particle].push_back(l);
   }


   std::vector<std::string>::const_iterator itStr;


   //pbarp reaction
   _pbarpReaction=std::shared_ptr<pbarpReaction>(new pbarpReaction(_prodChannelInfoList, id,_lmax));

   //preparations for prod key replacements
   std::vector<std::string> replProdKeyVec = _thePbarpParser->replaceProdKey();
   std::map<std::string, std::string> repProdKeyNames;

   for ( itStr = replProdKeyVec.begin(); itStr != replProdKeyVec.end(); ++itStr){
      std::stringstream stringStr;
      stringStr << (*itStr);
      std::string oldStr;
      stringStr >> oldStr;

      std::string newStr;
      stringStr >> newStr;
      repProdKeyNames[oldStr]=newStr;
   }

   //fill prodDecayList
   std::vector<std::string> additionalStringVecDummy;
   std::string dynTypeDefault="WoDynamics";

   additionalStringVecDummy.push_back("0.197");
   if(_thePbarpParser->productionFormalism()=="Cano"){
      std::vector< std::shared_ptr<IsobarLSDecay> > prodDecs= _pbarpReaction->productionDecays();
      std::vector< std::shared_ptr<IsobarLSDecay> >::iterator itDec;
      for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
	//	if(_theParser->useProductionBarrier()) (*itDec)->enableProdBarrier(_theParser->qRProduction());
	//	else (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
	 _prodDecList->addDecay(*itDec);
      }
   }
   else if(_thePbarpParser->productionFormalism()=="Tensor"){
      std::vector< std::shared_ptr<IsobarTensorDecay> > prodDecs= _pbarpReaction->productionTensorDecays();
      std::vector< std::shared_ptr<IsobarTensorDecay> >::iterator itDec;
      for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
         _prodDecList->addDecay(*itDec);
	 //	 (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
      }
   }
   else if(_thePbarpParser->productionFormalism()=="Heli"){
      std::vector< std::shared_ptr<IsobarHeliDecay> > prodDecs= _pbarpReaction->productionHeliDecays();
      std::vector< std::shared_ptr<IsobarHeliDecay> >::iterator itDec;
      for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
	//	if(_theParser->useProductionBarrier()) (*itDec)->enableProdBarrier(_theParser->qRProduction());
	//	else (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
	 _prodDecList->addDecay(*itDec);
      }
   }
   else{
      Alert <<"production formalism\t" << _thePbarpParser->productionFormalism() << "\t is not supported!!!" << endmsg;
      exit(0);
   }

   //enable dynamics and replace mass keys
   std::map<std::string, std::string>::iterator itMapStrStr;
   for (itMapStrStr=repProdKeyNames.begin(); itMapStrStr!=repProdKeyNames.end(); ++itMapStrStr){
      _prodDecList->replaceProdKey(itMapStrStr->first, itMapStrStr->second);
   }

   std::vector< std::shared_ptr<AbsDecay> > theProdDecs=_prodDecList->getList();
   std::vector< std::shared_ptr<AbsDecay> >::iterator itAbsDec;
   for (itAbsDec=theProdDecs.begin(); itAbsDec!=theProdDecs.end(); ++itAbsDec){
     if((*itAbsDec)->prodChannelInfo()->withProdBarrier()) (*itAbsDec)->enableProdBarrier();
     else (*itAbsDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
   }


   //set prefactor for production and decay amplitudes
   AbsChannelEnv::setPrefactors();

  //replace suffixes for fit parameter
  AbsChannelEnv::replaceParameterSuffixes();

   //replace mass key for decays
  AbsChannelEnv::replaceMassKeys();

   //add dynamics
  AbsChannelEnv::addDynamics();

   //set decay levels
  AbsChannelEnv::setDecayLevels();

   // spin density particles
   _spinDensity = _thePbarpParser->spinDensityNames();

}


