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
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"

epemReaction::epemReaction(std::vector<std::pair<Particle*, Particle*> >& prodPairs, ChannelID channelID) :
   _channelID(channelID)
  ,_epemIGJPC(new IGJPC(1,-1,-1, 0, -1))
{
    std::vector<std::pair<Particle*, Particle*> >::iterator itPartPairs;
    for (itPartPairs=prodPairs.begin(); itPartPairs!= prodPairs.end(); ++itPartPairs){
      //      std::string decName=(*itJPC)->name();
      std::shared_ptr<IsobarLSDecay> currentDec(new IsobarLSDecay( _epemIGJPC, itPartPairs->first, itPartPairs->second, _channelID, "epem"));
      currentDec->extractStates();

      if (currentDec->JPCLSAmps().size()>0){
	_prodCanoDecs.push_back(currentDec);
	std::shared_ptr<IsobarHeliDecay> currentHeliDec(new IsobarHeliDecay( _epemIGJPC,itPartPairs->first, itPartPairs->second, _channelID, "epem"));
	currentHeliDec->extractStates();
	_prodHeliDecs.push_back(currentHeliDec);
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
