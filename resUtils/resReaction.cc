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

// resReaction class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>

#include "resUtils/resReaction.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/ProdChannelInfo.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"

resReaction::resReaction(Particle* motherParticle, std::vector<std::shared_ptr<ProdChannelInfo> > prodChannelInfoList, ChannelID channelID) :
   _channelID(channelID)
  ,_motherParticle(motherParticle)
{

  std::vector<std::shared_ptr<ProdChannelInfo> >::iterator itProd;
  for (itProd=prodChannelInfoList.begin(); itProd!= prodChannelInfoList.end(); ++itProd){
   std::pair<Particle*, Particle*> particlePair = (*itProd)->productionPair();

   std::shared_ptr<IsobarLSDecay> currentDec(new IsobarLSDecay( _motherParticle, particlePair.first, particlePair.second, _channelID));
      currentDec->setProductionAmp();
      currentDec->extractStates();

      if (currentDec->JPCLSAmps().size()>0){
	_prodCanoDecs.push_back(currentDec);
	std::shared_ptr<IsobarHeliDecay> currentHeliDec(new IsobarHeliDecay( _motherParticle, particlePair.first, particlePair.second, _channelID));
	currentHeliDec->setProductionAmp();
	currentHeliDec->extractStates();
	_prodHeliDecs.push_back(currentHeliDec);
      }
    }

}

resReaction::~resReaction(){
}

void resReaction::print(std::ostream& os) const{
  os << "\n res reaction\n";

  os << "\n ***** decay chains *******\n";
  std::vector< std::shared_ptr<IsobarLSDecay> >::const_iterator itIso;
  for( itIso=_prodCanoDecs.begin(); itIso!=_prodCanoDecs.end(); ++itIso){
    (*itIso)->print(os);
  }
}
