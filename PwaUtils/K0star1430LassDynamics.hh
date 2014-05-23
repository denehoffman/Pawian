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

// K0star1430LassDynamics class definition file. -*- C++ -*-
// Copyright 2014 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <memory>

#include "PwaUtils/AbsDynamics.hh"

class K0star1430Lass;

class K0star1430LassDynamics : public AbsDynamics{

public:
  K0star1430LassDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother);
  virtual ~K0star1430LassDynamics();

  virtual complex<double> eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom=0);
  
  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr);
  virtual bool checkRecalculation(fitParams& theParamVal);
  virtual void updateFitParams(fitParams& theParamVal);
  virtual void setMassKey(std::string& theMassKey);

protected:
  //  std::string _key;
  std::string _aLASSKey;
  std::string _rLASSKey;
  std::string _BLASSKey;
  std::string _phiBKey;
  std::string _RLASSKey;
  std::string _phiRKey;

  double _currentMass;
  double _currentWidth;
  double _currentaLASS;
  double _currentrLASS;
  double _currentBLASS;
  double _currentphiB;
  double _currentRLASS;
  double _currentphiR;

  std::map<int, complex<double> >  _cachedMap; 

private:

};
