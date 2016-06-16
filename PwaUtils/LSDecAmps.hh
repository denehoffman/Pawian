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

// LSDecAmps class definition file. -*- C++ -*-
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

class IsobarLSDecay;
class AbsDecay;
class AbsPawianParameters;

class LSDecAmps : public AbsXdecAmp{

public:

  // create/copy/destroy:

  ///Constructor
  LSDecAmps(std::shared_ptr<IsobarLSDecay> theDec, ChannelID channelID);
  LSDecAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID);
  /** Destructor */
  virtual ~LSDecAmps();


  // Getters:

  virtual void print(std::ostream& os) const;
  std::vector< std::shared_ptr<const LScomb> >& lsVec() {return _LSs;}

  // virtual void getDefaultParams(fitParCol& fitVal, fitParCol& fitErr);
  virtual complex<double> XdecAmp(Spin& lamX, EvtData* theData, AbsXdecAmp* grandmaAmp);
  virtual complex<double> XdecPartAmp(Spin& lamX, Spin& lamDec, short fixDaughterNr,
				      EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp);
  
  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);
  // virtual bool checkRecalculation(fitParCol& theParamVal);
  // virtual void updateFitParams(fitParCol& theParamVal);
  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void fillParamNameList();
  virtual void calcDynamics(EvtData* theData, AbsXdecAmp* grandmaAmp=0);
  // virtual void retrieveWignerDs(EvtData* theData);

protected:
  //  std::vector< std::shared_ptr<const JPCLS> > _JPCLSs;
  std::vector< std::shared_ptr<const LScomb> > _LSs;
  double _factorMag;
  double _parityFactor;
  std::map< std::shared_ptr<const LScomb>, map<Spin,map<Spin, double > >, pawian::Collection::SharedPtrLess > _cgPreFactor;
  std::map< Spin, map<Spin, map<std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > > > _cgPreFactor_LamLamLSMap;

  std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > _currentParamMags;
  std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > _currentParamPhis;
  std::map< std::shared_ptr<const LScomb>, complex<double>, pawian::Collection::SharedPtrLess > _currentParamMagExpi;
  std::map< std::shared_ptr<const LScomb>, complex<double>, pawian::Collection::SharedPtrLess > _currentParamPreFacMagExpi;

  Spin _Smax;
  //  std::map<std::thread::id, std::map<Spin, complex<double> > > _cachedDynLSMap;
  std::map<std::thread::id, std::map<unsigned short, std::map<Spin, complex<double> > > > _cachedDynIdLSMap;
  //  std::map<Id3StringType, complex<double> >* _currentWignerDMap;
  // std::map<Id3StringType, complex<double> > _currentWignerDMap;

  void  fillCgPreFactor();
  virtual complex<double> lsLoop(AbsXdecAmp* grandmaAmp, Spin& lamX, EvtData* theData, Spin& lam1Min, Spin& lam1Max, Spin& lam2Min, Spin& lam2Max, bool withDecs);


private:




};



