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

// GammapChannelEnv class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf, Julian Pychy

#pragma once

#include "PwaUtils/AbsChannelEnv.hh"


class gammapParser;
class gammapReaction;


class GammapChannelEnv : public AbsChannelEnv
{
public:
   virtual void setupChannel(ChannelID id);
   GammapChannelEnv(gammapParser* theParser);
  //  const float gammaMomentum() const {return _gammaMomentum;}
   const int lmax() const {return _lmax;}
   const double beamPolFraction() const {return _beamPolFraction;}
   const double beamPolAngle() const {return _beamPolAngle;}

  std::shared_ptr<gammapReaction> reaction() {return _gammapReaction;}
   std::vector<std::string>& spinDensityNames(){ return _spinDensity;}
   std::map<std::string, std::vector<short> > dropGammapLForParticleData(){ return _dropGammapLForParticleData; }
  virtual const std::string  channelTypeName() {return "gammap";}
  virtual std::shared_ptr<AbsHist> CreateHistInstance(std::string additionalSuffix);

private:
   int _lmax;
   double _beamPolFraction;
   double _beamPolAngle;
  //  float _gammaMomentum;
   gammapParser* _theGamPParser;

  std::shared_ptr<gammapReaction> _gammapReaction;
  std::vector<std::string> _spinDensity;
  std::map<std::string, std::vector<short> > _dropGammapLForParticleData;
};
