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

// resBaseLh class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <memory>
#include <boost/function.hpp>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "Minuit2/MnUserParameters.h"

class AbsXdecAmp;
class resReaction;
class LSDecAmps;

class resBaseLh : public AbsLh {

public:
  // resBaseLh(std::shared_ptr<const EvtDataBaseList>);
  resBaseLh(ChannelID channelID);

  virtual ~resBaseLh();

  virtual AbsLh* clone_() const{
    AbsLh* theClone=new resBaseLh(_channelID);
    theClone->setDataVec(_evtDataVec);
    theClone->setMcVec(_evtMCVec);
    return theClone;
  }

  virtual double calcEvtIntensity( EvtData* theData, fitParCol& theParamVal);
  virtual complex<double> calcProdPartAmp(Spin lamX, Spin lamDec, std::string nameDec, EvtData* theData,
					  std::map <std::shared_ptr<const JPCLS>,
					  std::vector< std::shared_ptr<AbsXdecAmp> >,
					  pawian::Collection::SharedPtrLess > pbarpAmps);

  virtual complex<double> calcSpinDensity(Spin M1, Spin M2, std::string& nameDec, EvtData* theData);

  virtual void print(std::ostream& os) const;


protected:
  std::shared_ptr<resReaction> _resReactionPtr;

  int _highestJFsp;
  bool _isHighestJaPhoton;
  Spin _Jmother;
private:

  void initialize();
};
