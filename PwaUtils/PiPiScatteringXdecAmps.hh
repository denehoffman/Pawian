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

// PiPiScatteringXdecAmps class definition file. -*- C++ -*-
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

class AbsDecay;
class AbsPawianParameters;
class PiPiScatteringDecay;
class TMatrixDynamics;

class PiPiScatteringXdecAmps : public AbsXdecAmp{

public:

  // create/copy/destroy:

  ///Constructor
  PiPiScatteringXdecAmps(std::shared_ptr<PiPiScatteringDecay> theDec, ChannelID channelID);
  PiPiScatteringXdecAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID);

  /** Destructor */
  virtual ~PiPiScatteringXdecAmps();


  // Getters:

  virtual void print(std::ostream& os) const;

  // virtual void getDefaultParams(fitParCol& fitVal, fitParCol& fitErr);
  virtual complex<double> XdecPartAmp(const Spin& lamX, Spin& lamDec, short fixDaughterNr,
                                      EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp=0);
  virtual complex<double> XdecAmp(const Spin& lamX, EvtData* theData, AbsXdecAmp* grandmaAmp=0);
  
  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void fillParamNameList();
  virtual void calcDynamics(EvtData* theData, AbsXdecAmp* grandmaAmp=0);
  void doTcheck();
  
protected:

  virtual void initialize();
  std::shared_ptr<TMatrixDynamics> _tMatDyn;
  int _orbMom;
  
private:




};



