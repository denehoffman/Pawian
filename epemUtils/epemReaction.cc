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

// epemReaction class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>

#include "epemUtils/epemReaction.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/IsobarTensorDecay.hh"
#include "PwaUtils/IsobarTensorPsiToGamXDecay.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"

epemReaction::epemReaction(std::vector<std::pair<Particle*, Particle*> >& prodPairs, ChannelID channelID) :
   _channelID(channelID)
  ,_epemIGJPC(new IGJPC(1,-1,-1, 0, -1))
{
    std::vector<std::pair<Particle*, Particle*> >::iterator itPartPairs;
    for (itPartPairs=prodPairs.begin(); itPartPairs!= prodPairs.end(); ++itPartPairs){
      std::shared_ptr<IsobarLSDecay> currentDec(new IsobarLSDecay( _epemIGJPC, itPartPairs->first, itPartPairs->second, _channelID, "epem"));
      currentDec->setProductionAmp();
      currentDec->extractStates();

      if (currentDec->JPCLSAmps().size()>0){
	_prodCanoDecs.push_back(currentDec);
	std::shared_ptr<IsobarTensorDecay> currentTensorDec;

	currentTensorDec=std::shared_ptr<IsobarTensorDecay> (new IsobarTensorDecay( _epemIGJPC,itPartPairs->first, itPartPairs->second, _channelID, "epem")); //workaround
	
	currentTensorDec->setProductionAmp();
	currentTensorDec->extractStates();	
	_prodTensorDecs.push_back(currentTensorDec);

	//fill Zou amplitudes
	std::shared_ptr<IsobarTensorDecay> currentTensorZouDec;
	if(itPartPairs->first->name()=="photon" || itPartPairs->second->name()=="photon"){
	  if(itPartPairs->first->name()=="photon"){
	    currentTensorZouDec=std::shared_ptr<IsobarTensorDecay> (new IsobarTensorPsiToGamXDecay( _epemIGJPC,itPartPairs->first, itPartPairs->second, _channelID, "epem"));
	  }
	  else{ //itPartPairs->second->name()=="photon"
	    currentTensorZouDec=std::shared_ptr<IsobarTensorPsiToGamXDecay> (new IsobarTensorPsiToGamXDecay( _epemIGJPC,itPartPairs->second, itPartPairs->first, _channelID, "epem"));
	  }
	  currentTensorZouDec->setProductionAmp();
	  currentTensorZouDec->extractStates();	
	  _prodTensorZouDecs.push_back(currentTensorZouDec);
	}
	else _prodTensorZouDecs.push_back(currentTensorDec); 

	std::shared_ptr<IsobarHeliDecay> currentHeliDec(new IsobarHeliDecay( _epemIGJPC,itPartPairs->first, itPartPairs->second, _channelID, "epem"));
	currentHeliDec->setProductionAmp();
	currentHeliDec->extractStates();
	_prodHeliDecs.push_back(currentHeliDec);

	if(itPartPairs->first->name()=="photon" || itPartPairs->second->name()=="photon"){
	  std::shared_ptr<IsobarHeliDecay> currentHeliMultipoleDec;
	  if(itPartPairs->first->name()=="photon"){
	    currentHeliMultipoleDec=std::shared_ptr<IsobarHeliDecay>(new IsobarHeliDecay( _epemIGJPC,itPartPairs->first, itPartPairs->second, _channelID, "epem", "IsobarHeliMultipoleDecay"));
	  }
	  else{
	    currentHeliMultipoleDec=std::shared_ptr<IsobarHeliDecay>(new IsobarHeliDecay( _epemIGJPC,itPartPairs->second, itPartPairs->first, _channelID, "epem", "IsobarHeliMultipoleDecay"));
	  }
	    currentHeliMultipoleDec->setProductionAmp();
	    currentHeliMultipoleDec->extractStates();
	    _prodHeliMultipoleDecs.push_back(currentHeliMultipoleDec); 
	}
	else _prodHeliMultipoleDecs.push_back(currentHeliDec); 
      }
    }

}

epemReaction::~epemReaction(){
}

void epemReaction::print(std::ostream& os) const{
  os << "\n e+e- reaction\n";

  os << "\n ***** decay chains *******\n";
  std::vector< std::shared_ptr<IsobarLSDecay> >::const_iterator itIso;
  for( itIso=_prodCanoDecs.begin(); itIso!=_prodCanoDecs.end(); ++itIso){
    (*itIso)->print(os);
  }
}
