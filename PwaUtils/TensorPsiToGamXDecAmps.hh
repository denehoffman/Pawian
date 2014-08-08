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

// TensorPsiToGamXDecAmps class definition file. -*- C++ -*-
// Copyright 2014 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <memory>

#include "PwaUtils/AbsXdecAmp.hh"

class IsobarTensorPsiToGamXDecay;
class AbsDecay;

class TensorPsiToGamXDecAmps : public AbsXdecAmp{

public:

  // create/copy/destroy:

  ///Constructor
  TensorPsiToGamXDecAmps(std::shared_ptr<IsobarTensorPsiToGamXDecay> theDec, ChannelID channelID);
  TensorPsiToGamXDecAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID);
  /** Destructor */
  virtual ~TensorPsiToGamXDecAmps();


  // Getters:

  virtual complex<double> XdecAmp(Spin& lamX, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp);
  virtual complex<double> XdecPartAmp(Spin& lamX, Spin& lamDec, short fixDaughterNr,
				      EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp);

  virtual void print(std::ostream& os) const;

  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr);
  virtual bool checkRecalculation(fitParams& theParamVal);
  virtual void updateFitParams(fitParams& theParamVal);

protected:
  std::map< int, double> _currentParamMags;
  std::map< int, double> _currentParamPhis;
  short _noOfAmps;
  std::vector<std::string> _MagParamNames;
  std::vector<std::string> _PhiParamNames;
  std::vector<double> _currentParamLocalMags;
  std::vector<double> _currentParamLocalPhis;

  complex<double> daughterAmp(Spin& lam2, EvtData* theData, Spin& lamFs);

private:




};



