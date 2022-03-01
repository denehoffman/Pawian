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

// LSOmegaTo3PiDecAmps class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <memory>

//#include "PwaUtils/LSDecRefAmps.hh"
#include "PwaUtils/LSDecAmps.hh"
#include "PwaUtils/OmegaTo3PiLSDecay.hh"
//class OmegaTo3PiLSDecay;

//class LSOmegaTo3PiDecAmps : public LSDecRefAmps{
class LSOmegaTo3PiDecAmps : public LSDecAmps{
public:

  // create/copy/destroy:

  ///Constructor
  LSOmegaTo3PiDecAmps(std::shared_ptr<OmegaTo3PiLSDecay> theDec, ChannelID channelID);

  /** Destructor */
  virtual ~LSOmegaTo3PiDecAmps();


  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar);
  // Getters:
  
  virtual complex<double> XdecAmp(const Spin& lamX, EvtData* theData, AbsXdecAmp* grandmaAmp);
  virtual complex<double> XdecPartAmp(const Spin& lamX, Spin& lamDec, short fixDaughterNr,
				      EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp);

  virtual void print(std::ostream& os) const;

protected:

  std::string _lambdaDecKey;
  std::map< std::shared_ptr<const LScomb>, complex<double>, pawian::Collection::SharedPtrLess > _currentParamMagExpi;
private:




};



