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

// BreitWignerDynamics class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <memory>

#include "PwaUtils/BreitWignerDynamics.hh"

class BreitWignerRelDynamics : public BreitWignerDynamics{

public:
  BreitWignerRelDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, double massSumDaughter1, double massSumDaughter2);
  virtual ~BreitWignerRelDynamics();

  virtual std::string type() {return "BreitWignerRelDynamics";}
  virtual complex<double> eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom=0);
  

protected:
  double _fsp1Mass;
  double _fsp2Mass;
private:

};
