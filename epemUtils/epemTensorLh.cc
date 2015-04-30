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

// epemTensorLh class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "epemUtils/epemBaseLh.hh"
#include "epemUtils/epemTensorLh.hh"
#include "epemUtils/epemReaction.hh"
#include "epemUtils/EpemChannelEnv.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarTensorDecay.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


epemTensorLh::epemTensorLh(ChannelID channelID) : 
  epemBaseLh(channelID)
{
  initialize();
}


epemTensorLh::~epemTensorLh()
{;
}


void epemTensorLh::print(std::ostream& os) const{

}


void  epemTensorLh::initialize(){
  epemBaseLh::initialize();
  //  std::vector< std::shared_ptr<IsobarTensorDecay> > theDecs = _epemReactionPtr->productionTensorDecays();
  // std::vector< std::shared_ptr<IsobarTensorDecay> >::iterator it;
  // for (it=theDecs.begin(); it!=theDecs.end(); ++it){
  //   std::shared_ptr<AbsXdecAmp> currentAmp=XdecAmpRegistry::instance()->getXdecAmp(_channelID, (*it)->absDecPtr());
  //   _decAmps.push_back(currentAmp);
  // }

  // std::vector< std::shared_ptr<AbsDecay> > theDecs = _epemChannelEnv->prodDecayList()->getList();
  // std::vector< std::shared_ptr<AbsDecay> >::iterator it;
  // for (it=theDecs.begin(); it!=theDecs.end(); ++it){
  //   std::shared_ptr<AbsXdecAmp> currentAmp=XdecAmpRegistry::instance()->getXdecAmp(_channelID, (*it)->absDecPtr());
  //   _decAmps.push_back(currentAmp);
  // }
}




