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

// IsobarTensorDecay class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <memory>

#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/DataUtils.hh"
#include "Utils/PawianCollectionUtils.hh"

class Particle;
class EvtData;

class IsobarTensorDecay : public IsobarLSDecay{

public:
  IsobarTensorDecay(Particle* mother, Particle* daughter1, Particle* daughter2, ChannelID channelID);
  IsobarTensorDecay(std::shared_ptr<const IGJPC> motherIGJPCPtr, Particle* daughter1, Particle* daughter2, ChannelID channelID, std::string motherName="pbarp");
  virtual ~IsobarTensorDecay();

  virtual void fillWignerDs(std::map<std::string , Vector4<double> >& fsMap, Vector4<double>& prodParticle4Vec, EvtData* evtData);
  virtual void print(std::ostream& os) const;
  virtual std::string type() {return "IsobarTensorDecay";}

protected:
  PolVector _polMother;
  PolVector _polDaughter1;
  PolVector _polDaughter2;
  LeviCivitaTensor _lctTensor;
  MetricTensor _metricTensor;
};
