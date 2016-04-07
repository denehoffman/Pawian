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

// TensorOmegaTo3PiDecAmps class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <memory>

#include "PwaUtils/AbsXdecAmp.hh"

class OmegaTo3PiTensorDecay;
class Particle;
class AbsPawianParameters;

typedef CacheVector<std::map<unsigned short,std::map<Id1StringType, complex<float> > > >  tensorOmegaTo3PiCachedIdMap;

class TensorOmegaTo3PiDecAmps : public AbsXdecAmp{

public:

  // create/copy/destroy:

  ///Constructor
  TensorOmegaTo3PiDecAmps(std::shared_ptr<OmegaTo3PiTensorDecay> theDec, ChannelID channelID);

  /** Destructor */
  virtual ~TensorOmegaTo3PiDecAmps();


  // Getters:

  virtual complex<double> XdecAmp(Spin& lamX, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp);
  virtual complex<double> XdecPartAmp(Spin& lamX, Spin& lamDec, short fixDaughterNr,
				      EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp);

  virtual void print(std::ostream& os) const;
  std::vector< std::shared_ptr<const LScomb> >& lsVec() {return _LSs;}

  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void fillParamNameList();

protected:
  std::vector< std::shared_ptr<const LScomb> > _LSs;
  double _factorMag;
  // std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > _currentParamMags;
  // std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > _currentParamPhis;
  std::map< std::shared_ptr<const LScomb>, complex<double>, pawian::Collection::SharedPtrLess > _currentMagExpi;

  Particle* _daughter1;
  Particle* _daughter2;
  Particle* _daughter3;
  tensorOmegaTo3PiCachedIdMap _cachedLocalAmpIdMap;
private:




};



