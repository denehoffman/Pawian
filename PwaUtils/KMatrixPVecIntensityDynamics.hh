//************************************************************************//
//									  //
//  Copyright 2019 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// KMatrixPVecIntensityDynamics class definition file. -*- C++ -*-
// Copyright 2019 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <memory>

#include "PwaUtils/KMatrixDynamics.hh"

class KMatrixPVecIntensityDynamics : public KMatrixDynamics{

public:
  KMatrixPVecIntensityDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::string& pathToConfigParser,  ChannelID channelID, std::string projectionParticleNames="");
  virtual ~KMatrixPVecIntensityDynamics();

  virtual std::string type() {return "KMatrixPVecIntensityDynamics";}
  virtual complex<double> eval(EvtData* theData, AbsXdecAmp* grandmaAmp=0, Spin OrbMom=0);
  virtual bool checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld); 
  virtual void addGrandMa(std::shared_ptr<AbsDecay> theDec) { return;} //dummy
protected:
  int _projectionCompareIndex;
  std::string _nameOfFVector;

private:

};
