//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// ProdParamDynamics class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <memory>

#include "PwaUtils/AbsDynamics.hh"

class ProdParamDynamics : public AbsDynamics{

public:
  ProdParamDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, ChannelID channelID, std::string type);
  virtual ~ProdParamDynamics();

  virtual std::string type() {return "ProdParamDynamics";}
  virtual complex<double> eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom=0);
  
  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);
  // virtual bool checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld);
  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void fillParamNameList();
  virtual void fillMasses(EvtData* theData);

protected:

  vector<double> _currentPolParams;
  double _currentExpParam;
  vector<std::string> _fitPolParNames;
  std::string _fitExpParName;
  unsigned int _polOrder;
  ChannelID _channelID;

private:

};
