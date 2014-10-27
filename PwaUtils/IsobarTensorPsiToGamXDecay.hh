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

// IsobarTensorPsiToGamXDecay class definition file. -*- C++ -*-
// Copyright 2014 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <memory>

#include "PwaUtils/IsobarTensorDecay.hh"

class Particle;
class EvtData;

class IsobarTensorPsiToGamXDecay : public IsobarTensorDecay{

public:
  IsobarTensorPsiToGamXDecay(Particle* mother, Particle* daughter1, Particle* daughter2, ChannelID channelID);
  IsobarTensorPsiToGamXDecay(std::shared_ptr<const IGJPC> motherIGJPCPtr, Particle* daughter1, Particle* daughter2, ChannelID channelID, std::string motherName="pbarp");
  virtual ~IsobarTensorPsiToGamXDecay();

  virtual void fillWignerDs(std::map<std::string , Vector4<double> >& fsMap, Vector4<double>& prodParticle4Vec, EvtData* evtData, std::string& refKey);
  virtual void print(std::ostream& os) const;
  virtual std::string type() {return "IsobarTensorPsiToGamXDecay";}

  virtual void enableProdBarrier(double qRValue);
  virtual short noOfAmplitudes() const {return _noOfAmps;}
  virtual std::map<int, int> ampLMap() {return _ampLMap;}
protected:
  bool _XisEven;
  short _noOfAmps;
  std::map<int, int> _ampLMap;

  void fillAmpLMap();
};
