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

#pragma once

#include "PwaUtils/AbsChannelEnv.hh"


class pbarpParser;
class pbarpReaction;


class PbarpChannelEnv : public AbsChannelEnv
{
public:
  void setup(ChannelID id);
  PbarpChannelEnv(pbarpParser* theParser);
  const float pbarMomentum() const {return _pbarMomentum;}
  const int lmax() const {return _lmax;}
  
  std::shared_ptr<pbarpReaction> reaction() {return _pbarpReaction;}
  std::vector<std::string>& spinDensityNames(){ return _spinDensity;}
  std::map<std::string, std::vector<short> > dropPbarpLForParticleData(){ return _dropPbarpLForParticleData; }
  virtual const std::string  channelTypeName() {return "pbarp";}
  virtual std::shared_ptr<AbsHist> CreateHistInstance(std::string additionalSuffix);
  
private:
  int _lmax;
  float _pbarMomentum;
  pbarpParser* _thePbarpParser;
  
  std::shared_ptr<pbarpReaction> _pbarpReaction;
  std::vector<std::string> _spinDensity;
  std::map<std::string, std::vector<short> > _dropPbarpLForParticleData;
};
