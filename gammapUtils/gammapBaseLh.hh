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

// gammapBaseLh class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

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

class AbsXdecAmp;
class gammapReaction;
class LSDecAmps;
class AbsPawianParameters;

class gammapBaseLh : public AbsLh {

public:
  // gammapBaseLh(std::shared_ptr<const EvtDataBaseList>);
  gammapBaseLh(ChannelID channelID);

  virtual ~gammapBaseLh();

  virtual AbsLh* clone_() const{
    AbsLh* theClone=new gammapBaseLh(_channelID);
    theClone->setDataVec(_evtDataVec);
    theClone->setMcVec(_evtMCVec);
    return theClone;
  }

  // virtual double calcEvtIntensity( EvtData* theData, fitParCol& theParamVal);
  virtual complex<double> calcProdPartAmp(Spin lamX, Spin lamDec, std::string nameDec, EvtData* theData,
					  std::map <std::shared_ptr<const JPCLS>,
					  std::vector< std::shared_ptr<AbsXdecAmp> >,
					  pawian::Collection::SharedPtrLess > gammapAmps);

  virtual complex<double> calcSpinDensity(Spin M1, Spin M2, std::string& nameDec, EvtData* theData);

  // virtual void getDefaultParams(fitParCol& fitVal, fitParCol& fitErr);
  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);
  
  // virtual void updateFitParams(fitParCol& theParamVal);
  
  virtual void print(std::ostream& os) const;


protected:
  std::shared_ptr<gammapReaction> _gammapReactionPtr;
  std::map< std::shared_ptr<const jpcRes>, std::vector< std::shared_ptr<const IGJPC> >, pawian::Collection::SharedPtrLess> _jpcToIGJPCMap;
  std::map< std::shared_ptr<const jpcRes>, std::vector< std::shared_ptr<const JPCLS> >, pawian::Collection::SharedPtrLess> _jpcToJPCljMap;

  std::map< std::shared_ptr<const jpcRes>, double, pawian::Collection::SharedPtrLess > _currentParamMags;
  std::map< std::shared_ptr<const jpcRes>, double, pawian::Collection::SharedPtrLess > _currentParamPhis;

  std::map< std::shared_ptr<const jpcRes>, double, pawian::Collection::SharedPtrLess> _currentParamJPCIsos12;
  std::map< std::shared_ptr<const jpcRes>, double, pawian::Collection::SharedPtrLess> _currentParamJPCIsos32;
  std::vector< std::shared_ptr<AbsXdecAmp> > _iso12DecAmps;
  std::vector< std::shared_ptr<AbsXdecAmp> > _iso32DecAmps;

  virtual void fillIsos();

private:

  void initialize();
};
