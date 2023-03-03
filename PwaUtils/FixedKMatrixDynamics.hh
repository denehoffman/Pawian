//************************************************************************//
//									  //
//  Copyright 2023 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// FixedKMatrixDynamics class definition file. -*- C++ -*-
// Copyright 2023 Bertram Kopf

#pragma once

#include "PwaUtils/KMatrixDynamics.hh"

class AbsPawianParameters;

class FixedKMatrixDynamics : public KMatrixDynamics{

public:
  FixedKMatrixDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::string& pathToConfigParser, ChannelID channelID, std::string projectionParticleNames="");

  virtual ~FixedKMatrixDynamics();

  virtual std::string type() {return "FixedKMatrixDynamics";}
  
  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void fillParamNameList();
  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar);
 
protected:
  std::shared_ptr<AbsPawianParameters> _kMatFixParams;

private:

};
