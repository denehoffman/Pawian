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

// FormXDecAmp class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <memory>

#include "PwaUtils/AbsXdecAmp.hh"

class FormationDecay;
class AbsDecay;
class AbsPawianParameters;

class FormXDecAmps : public AbsXdecAmp{

public:

  // create/copy/destroy:

  ///Constructor
  FormXDecAmps(std::shared_ptr<FormationDecay> theDec, ChannelID channelID);
  FormXDecAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID);
  /** Destructor */
  virtual ~FormXDecAmps();


  // Getters:

  virtual void print(std::ostream& os) const;

  virtual complex<double> XdecAmp(Spin& lamX, EvtData* theData, AbsXdecAmp* grandmaAmp);
  virtual complex<double> XdecPartAmp(Spin& lamX, Spin& lamDec, short fixDaughterNr,
				      EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp);
  virtual void fillParamNameList();
  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void calcDynamics(EvtData* theData, AbsXdecAmp* grandmaAmp=0);

  virtual complex<double> daughterAmp(Spin& lam1, Spin& lam2, EvtData* theData);
  virtual void cacheAmplitudes();
  virtual bool checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld);

protected:
  std::map<Spin, std::string> _paramNameMap;
  std::map<Spin, complex<double> > _currentParamMap;

  virtual void initialize();

private:




};



