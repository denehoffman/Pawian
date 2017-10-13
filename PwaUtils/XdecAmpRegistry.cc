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

// XdecAmpRegistry class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/XdecAmpRegistry.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/HeliMultipoleDecNonRefAmps.hh"
#include "PwaUtils/IsobarTensorDecay.hh"
#include "PwaUtils/IsobarTensorPsiToGamXDecay.hh"
#include "PwaUtils/OmegaTo3PiLSDecay.hh"
#include "PwaUtils/OmegaTo3PiTensorDecay.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/LSDecAmps.hh"
#include "PwaUtils/HeliDecAmps.hh"
#include "PwaUtils/TensorDecAmps.hh"
#include "PwaUtils/TensorPsiToGamXDecAmps.hh"
#include "PwaUtils/LSOmegaTo3PiDecAmps.hh"
#include "PwaUtils/TensorOmegaTo3PiDecAmps.hh"
#include "PwaUtils/PiPiScatteringDecay.hh"
#include "PwaUtils/PiPiScatteringXdecAmps.hh"
#include "PwaUtils/FormationDecay.hh"
#include "PwaUtils/FormXDecAmps.hh"
#include "ErrLogger/ErrLogger.hh"


XdecAmpRegistry* XdecAmpRegistry::_instance=0;

XdecAmpRegistry* XdecAmpRegistry::instance()
{
  if (0==_instance) _instance = new XdecAmpRegistry();
  return _instance;
}

XdecAmpRegistry::XdecAmpRegistry()
{
}

XdecAmpRegistry::~XdecAmpRegistry()
{
}

std::shared_ptr<AbsXdecAmp> XdecAmpRegistry::getXdecAmp(short channelID, std::shared_ptr<AbsDecay> theAbsDecay){

  std::shared_ptr<AbsXdecAmp> result;

  std::string theName=theAbsDecay->name();
  std::map<std::string, std::shared_ptr<AbsXdecAmp> >::iterator it = _xDecAmpMap[channelID].find(theName);
  if (it !=_xDecAmpMap[channelID].end()) result=it->second;
  else{
    if(theAbsDecay->type()=="IsobarLSDecay"){
      std::shared_ptr<IsobarLSDecay> decLS =  std::dynamic_pointer_cast<IsobarLSDecay>(theAbsDecay);
      result=std::shared_ptr<AbsXdecAmp>(new LSDecAmps(decLS, channelID));
    }
    else if(theAbsDecay->type()=="IsobarHeliDecay"){
      std::shared_ptr<IsobarHeliDecay> decLamLam =  std::dynamic_pointer_cast<IsobarHeliDecay>(theAbsDecay);
      result=std::shared_ptr<AbsXdecAmp>(new HeliDecAmps(decLamLam, channelID));
    }
 else if(theAbsDecay->type()=="IsobarHeliMultipoleDecay"){
   std::shared_ptr<IsobarHeliDecay> decLamLam =  std::dynamic_pointer_cast<IsobarHeliDecay>(theAbsDecay);
   if(theAbsDecay->whichDecayLevel()==AbsDecay::decayLevel::isProdAmp || theAbsDecay->whichDecayLevel()==AbsDecay::decayLevel::firstLevel) result=std::shared_ptr<AbsXdecAmp>(new HeliMultipoleDecNonRefAmps(decLamLam, channelID));
   else{
     Alert << "creation HeliMultipoleDecNonRefAmps not supported for non-production amps!!!" << endmsg;
     exit(0); 
   }
  }
    else if(theAbsDecay->type()=="IsobarTensorDecay"){
      std::shared_ptr<IsobarTensorDecay> decTensor =  std::dynamic_pointer_cast<IsobarTensorDecay>(theAbsDecay);
      result=std::shared_ptr<AbsXdecAmp>(new TensorDecAmps(decTensor, channelID));
    }
    else if(theAbsDecay->type()=="IsobarTensorPsiToGamXDecay"){
      std::shared_ptr<IsobarTensorPsiToGamXDecay> decTensor =  std::dynamic_pointer_cast<IsobarTensorPsiToGamXDecay>(theAbsDecay);
      result=std::shared_ptr<AbsXdecAmp>(new TensorPsiToGamXDecAmps(decTensor, channelID));
    }
    else if(theAbsDecay->type()=="OmegaTo3PiLSDecay"){
      std::shared_ptr<OmegaTo3PiLSDecay> decOmega =  std::dynamic_pointer_cast<OmegaTo3PiLSDecay>(theAbsDecay);
      result=std::shared_ptr<AbsXdecAmp>(new LSOmegaTo3PiDecAmps(decOmega, channelID));
    }
    else if(theAbsDecay->type()=="OmegaTo3PiTensorDecay"){
      std::shared_ptr<OmegaTo3PiTensorDecay> decOmega =  std::dynamic_pointer_cast<OmegaTo3PiTensorDecay>(theAbsDecay);
      result=std::shared_ptr<AbsXdecAmp>(new TensorOmegaTo3PiDecAmps(decOmega, channelID));
    }
    else if(theAbsDecay->type()=="PiPiScatteringDecay"){
      std::shared_ptr<PiPiScatteringDecay> pipiScatteringDecay =  std::dynamic_pointer_cast<PiPiScatteringDecay>(theAbsDecay);
      result=std::shared_ptr<AbsXdecAmp>(new PiPiScatteringXdecAmps(pipiScatteringDecay, channelID));
    }
    else if(theAbsDecay->type()=="FormationDecay"){
      std::shared_ptr<FormationDecay> formationDecay =  std::dynamic_pointer_cast<FormationDecay>(theAbsDecay);
      result=std::shared_ptr<FormXDecAmps>(new FormXDecAmps(formationDecay, channelID));
    }
    else{
      Alert << "can nor create XdecAmp object for theAbsDecay->name():\t" << theAbsDecay->name() << endmsg;
      exit(1);
    }
    result->fillParamNameList();
    _xDecAmpMap[channelID][result->name()]=result;
  }
  return result;
}

