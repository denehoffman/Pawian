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

// pbarpBaseLh class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/DataUtils.hh"
#include "Minuit2/MnUserParameters.h"

class AbsXdecAmp;
class pbarpReaction;
class LSDecAmps;

class pbarpBaseLh : public AbsLh {

public:
  // pbarpBaseLh(boost::shared_ptr<const EvtDataBaseList>);

  pbarpBaseLh();  

  virtual ~pbarpBaseLh();
  
  virtual AbsLh* clone_() const =0;

  virtual double calcEvtIntensity( EvtData* theData, fitParams& theParamVal);
  virtual complex<double> calcProdPartAmp(Spin lamX, Spin lamDec, std::string nameDec, EvtData* theData, 
					  std::map <boost::shared_ptr<const JPCLS>,
					  std::vector< boost::shared_ptr<AbsXdecAmp> >,
					  pawian::Collection::SharedPtrLess > pbarpAmps);

  virtual complex<double> calcSpinDensity(Spin M1, Spin M2, std::string& nameDec, EvtData* theData);
  virtual double calcSpinDensityNorm(std::string& nameDec, EvtData* theData);
  
  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr);
  virtual void updateFitParams(fitParams& theParamVal);
  
  virtual void print(std::ostream& os) const;
  
protected:
  boost::shared_ptr<pbarpReaction> _pbarpReactionPtr;
  std::vector< boost::shared_ptr<const JPCLS> > _jpclsStates;

  std::map <boost::shared_ptr<const JPCLS>, std::vector< boost::shared_ptr<AbsXdecAmp> >, pawian::Collection::SharedPtrLess > _decAmpsSinglet;
  std::map <boost::shared_ptr<const JPCLS>, std::vector< boost::shared_ptr<AbsXdecAmp> >, pawian::Collection::SharedPtrLess > _decAmpsTriplet0;
  std::map <boost::shared_ptr<const JPCLS>, std::vector< boost::shared_ptr<AbsXdecAmp> >, pawian::Collection::SharedPtrLess > _decAmpsTripletp1;
  std::map <boost::shared_ptr<const JPCLS>, std::vector< boost::shared_ptr<AbsXdecAmp> >, pawian::Collection::SharedPtrLess > _decAmpsTripletm1;
  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > _currentParamMags;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > _currentParamPhis;
  int _highestJFsp;
  bool _isHighestJaPhoton;
  
  virtual void fillMap(std::vector< boost::shared_ptr<const JPCLS> >& pbarpLSs, std::vector<boost::shared_ptr<AbsXdecAmp> >& decs, std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<AbsXdecAmp> >, pawian::Collection::SharedPtrLess >& toFill); 

private:

  void initialize();
};
