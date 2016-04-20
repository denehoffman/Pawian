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

// HeliDecAmps class definition file. -*- C++ -*-
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

class AbsPawianParameters;
class IsobarHeliDecay;
class AbsDecay;

class HeliDecAmps : public AbsXdecAmp{

public:

  // create/copy/destroy:

  ///Constructor
  HeliDecAmps(std::shared_ptr<IsobarHeliDecay> theDec, ChannelID channelID);
  HeliDecAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID);
  /** Destructor */
  virtual ~HeliDecAmps();


  // Getters:

  // virtual void getDefaultParams(fitParCol& fitVal, fitParCol& fitErr);
  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void print(std::ostream& os) const;
  // virtual bool checkRecalculation(fitParCol& theParamVal);
  // std::shared_ptr<const jpcRes>& jpcPtr() {return _JPCPtr;}
  std::vector< std::shared_ptr<const JPClamlam> >& jpclamlamVec() {return _JPClamlams;}
  // virtual void updateFitParams(fitParCol& theParamVal);
  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void fillParamNameList();

  virtual complex<double> XdecAmp(Spin& lamX, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp);
  virtual complex<double> XdecPartAmp(Spin& lamX, Spin& lamDec, short fixDaughterNr,
                                      EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp);
  virtual void calcDynamics(EvtData* theData, AbsXdecAmp* grandmaAmp);

protected:
  std::vector< std::shared_ptr<const JPClamlam> > _JPClamlams;
  double _factorMag;
  double _parityFactor;
  bool _isWeakDecay;
  std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > _currentParamMagLamLams;
  std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > _currentParamPhiLamLams;
  std::map< std::shared_ptr<const JPClamlam>, complex<double>, pawian::Collection::SharedPtrLess > _currentParamMagExpi;
  std::map< std::shared_ptr<const JPClamlam>, complex<double>, pawian::Collection::SharedPtrLess > _currentParamPreFacMagExpi;
  std::map< std::shared_ptr<const JPClamlam>, std::vector< std::shared_ptr<const JPClamlam> >, pawian::Collection::SharedPtrLess > _JPClamlamSymMap;
  std::map<std::thread::id, complex<double> > _cachedDynLMap;

  virtual void printCurrentAmpParams(Spin& lamX, Spin& lamFs);

private:




};



