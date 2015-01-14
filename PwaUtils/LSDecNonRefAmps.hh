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

// LSDecNonRefAmps class definition file. -*- C++ -*-
// Copyright 2014 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <memory>

#include "PwaUtils/AbsLSDecAmps.hh"

class IsobarLSDecay;
class AbsDecay;

class LSDecNonRefAmps : public AbsLSDecAmps{

public:

  // create/copy/destroy:

  ///Constructor
  LSDecNonRefAmps(std::shared_ptr<IsobarLSDecay> theDec, ChannelID channelID);
  LSDecNonRefAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID);
  /** Destructor */
  virtual ~LSDecNonRefAmps();


  // Getters:

  virtual complex<double> XdecAmp(Spin& lamX, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp);
  virtual complex<double> XdecPartAmp(Spin& lamX, Spin& lamDec, short fixDaughterNr,
				      EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp);

  virtual void calcDynamics(EvtData* theData, AbsXdecAmp* grandmaAmp=0);
protected:
  virtual complex<double> lsLoop(AbsXdecAmp* grandmaAmp, Spin& lamX, EvtData* theData, Spin& lam1Min, Spin& lam1Max, Spin& lam2Min, Spin& lam2Max, bool withDecs, Spin lamFs=0 );

  Spin _Smax;
  std::map<std::thread::id, std::map<Spin, complex<double> > > _cachedDynLSMap;

private:




};



