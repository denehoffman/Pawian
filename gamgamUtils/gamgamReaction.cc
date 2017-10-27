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

// gamgamReaction class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#include <getopt.h>
#include <fstream>

#include "gamgamUtils/gamgamReaction.hh"
#include "gamgamUtils/gamgamStates.hh"
#include "PwaUtils/FormationDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/ProdChannelInfo.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"

gamgamReaction::gamgamReaction(std::vector<std::shared_ptr<ProdChannelInfo> > prodChannelInfoList, ChannelID channelID, int jmax) :
   _channelID(channelID)
   ,_jmax(jmax)
   ,_gamgamStates(new gamgamStates(_jmax))
{
  std::vector<std::string> additionalStringVecDummy;  
     std::vector<std::shared_ptr<ProdChannelInfo> >::iterator itProd;
    for (itProd=prodChannelInfoList.begin(); itProd!= prodChannelInfoList.end(); ++itProd){
      if( (*itProd)->isFormation() ){
 	Particle* formationParticle=(*itProd)->formationParticle();
	if(0==formationParticle){
	  Alert << "formation particle does not exist!!!" << endmsg;
	  exit(1);
	}
	std::shared_ptr<const IGJPC> currentIGJPC=getIGJPCPtr(formationParticle);
	std::shared_ptr<FormationDecay> currentDec(new FormationDecay(currentIGJPC, formationParticle, _channelID, "gamgam"));
	currentDec->setProductionAmp();
	currentDec->setProdChannelInfo( *itProd );
	currentDec->extractStates();
	std::string currentFormDynType=(*itProd)->formationDynType();
	if(currentFormDynType !="non"){
	  if(currentFormDynType =="FormPol0" || currentFormDynType =="FormPol1" || currentFormDynType =="FormPol2"){
	  currentDec->enableDynamics(currentFormDynType, additionalStringVecDummy);
	  }
	}
	_prodFormationDecs.push_back(currentDec);
      }
      else{
	Alert << "gamma gamma production mode is not supported so far" << endmsg;
	exit(1);
      }
    }
}

gamgamReaction::~gamgamReaction(){
}

void gamgamReaction::print(std::ostream& os) const{
  os << "\n gamma gamma reaction\n";

  os << "\n ***** decay chains *******\n";
  // std::vector< std::shared_ptr<IsobarLSDecay> >::const_iterator itIso;
  // for( itIso=_prodCanoDecs.begin(); itIso!=_prodCanoDecs.end(); ++itIso){
  //   (*itIso)->print(os);
  // }
}
