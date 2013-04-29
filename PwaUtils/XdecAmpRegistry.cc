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
#include "PwaUtils/IsobarTensorDecay.hh"
#include "PwaUtils/OmegaTo3PiLSDecay.hh"
#include "PwaUtils/OmegaTo3PiTensorDecay.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/LSDecAmps.hh"
#include "PwaUtils/HeliDecAmps.hh"
#include "PwaUtils/TensorDecAmps.hh"
#include "PwaUtils/LSOmegaTo3PiDecAmps.hh"
#include "PwaUtils/TensorOmegaTo3PiDecAmps.hh"
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

boost::shared_ptr<AbsXdecAmp> XdecAmpRegistry::getXdecAmp(boost::shared_ptr<AbsDecay> theAbsXDec){

  boost::shared_ptr<AbsXdecAmp> result;

  std::string theName=theAbsXDec->name();
  std::map<std::string, boost::shared_ptr<AbsXdecAmp> >::iterator it = _xDecAmpMap.find(theName);
  if (it !=_xDecAmpMap.end()) result=it->second;
  else{
    if(theAbsXDec->type()=="IsobarLSDecay"){
      boost::shared_ptr<IsobarLSDecay> decLS =  boost::dynamic_pointer_cast<IsobarLSDecay>(theAbsXDec);
      result=boost::shared_ptr<AbsXdecAmp>(new LSDecAmps(decLS));
    }
    else if(theAbsXDec->type()=="IsobarHeliDecay"){
      boost::shared_ptr<IsobarHeliDecay> decLamLam =  boost::dynamic_pointer_cast<IsobarHeliDecay>(theAbsXDec);
      result=boost::shared_ptr<AbsXdecAmp>(new HeliDecAmps(decLamLam));
    }
    else if(theAbsXDec->type()=="IsobarTensorDecay"){
      boost::shared_ptr<IsobarTensorDecay> decTensor =  boost::dynamic_pointer_cast<IsobarTensorDecay>(theAbsXDec);
      result=boost::shared_ptr<AbsXdecAmp>(new TensorDecAmps(decTensor));
    }
    else if(theAbsXDec->type()=="OmegaTo3PiLSDecay"){
      boost::shared_ptr<OmegaTo3PiLSDecay> decOmega =  boost::dynamic_pointer_cast<OmegaTo3PiLSDecay>(theAbsXDec);
      result=boost::shared_ptr<AbsXdecAmp>(new LSOmegaTo3PiDecAmps(decOmega));
    }
    else if(theAbsXDec->type()=="OmegaTo3PiTensorDecay"){
      boost::shared_ptr<OmegaTo3PiTensorDecay> decOmega =  boost::dynamic_pointer_cast<OmegaTo3PiTensorDecay>(theAbsXDec);
      result=boost::shared_ptr<AbsXdecAmp>(new TensorOmegaTo3PiDecAmps(decOmega));
    }
    else{
      Alert << "can nor create XdecAmp object for theAbsXDec->name():\t" << theAbsXDec->name() << endmsg;
      exit(1); 
    }

    _xDecAmpMap[result->name()]=result;
  }
  return result;
}

