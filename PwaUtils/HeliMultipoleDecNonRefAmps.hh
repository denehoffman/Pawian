//************************************************************************//
//									  //
//  Copyright 2014 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// HeliMultipoleDecNonRefAmps class definition file. -*- C++ -*-
// Copyright 2015 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <memory>
#include <thread>

#include "PwaUtils/HeliDecNonRefAmps.hh"

class IsobarHeliDecay;
class AbsDecay;
class AbsPawianParameters;

class HeliMultipoleDecNonRefAmps : public HeliDecNonRefAmps{

public:

  // create/copy/destroy:

  ///Constructor
  HeliMultipoleDecNonRefAmps(std::shared_ptr<IsobarHeliDecay> theDec, ChannelID channelID);
  HeliMultipoleDecNonRefAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID);
  /** Destructor */
  virtual ~HeliMultipoleDecNonRefAmps();


  // Getters:
  virtual void getDefaultParams(fitParCol& fitVal, fitParCol& fitErr);
  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);
  //  virtual void print(std::ostream& os) const;
  virtual bool checkRecalculation(fitParCol& theParamVal);
  virtual void updateFitParams(fitParCol& theParamVal);

protected:

  std::map< int, double> _currentParamMags;
  std::map< int, double> _currentParamPhis;
  short _noOfAmps;
  std::vector<std::string> _MagParamNames;
  std::vector<std::string> _PhiParamNames;
  std::vector<double> _currentParamLocalMags;
  std::vector<double> _currentParamLocalPhis;
  std::vector< complex<double> > _currentParamLocalMagExpi;
  std::map<int, Spin> _JgammaMap;
  std::vector<Spin> _JgammaVec;
  Particle* _daughterGamma;
  Particle* _daughter2;
private:




};



