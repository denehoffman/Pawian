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

// pipiScatteringBaseLh class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

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
class PiPiScatteringXdecAmps;
class resReaction;


class pipiScatteringBaseLh : public AbsLh {

public:
  // pipiScatteringBaseLh(std::shared_ptr<const EvtDataBaseList>);
  pipiScatteringBaseLh(ChannelID channelID);

  virtual ~pipiScatteringBaseLh();

  virtual AbsLh* clone_() const{
    AbsLh* theClone=new pipiScatteringBaseLh(_channelID);
    theClone->setDataVec(_evtDataVec);
    theClone->setMcVec(_evtMCVec);
    return theClone;
  }

  virtual double calcLogLh(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual double addDataToLogLh(EvtData* dataEvt, std::shared_ptr<AbsPawianParameters> fitPar, LHData& theLHData);
  virtual double addMcToLogLh(EvtData* mcEvt, std::shared_ptr<AbsPawianParameters> fitPar, LHData& theLHData);
  virtual void calcLogLhDataClient(std::shared_ptr<AbsPawianParameters> fitPar, LHData& theLHData);
  virtual double calcEvtIntensity( EvtData* theData, std::shared_ptr<AbsPawianParameters> fitPar);
  virtual complex<double> calcProdPartAmp(Spin lamX, Spin lamDec, std::string nameDec, EvtData* theData,
					  std::map <std::shared_ptr<const JPCLS>,
					  std::vector< std::shared_ptr<AbsXdecAmp> >,
					  pawian::Collection::SharedPtrLess > resAmps);

  virtual complex<double> calcSpinDensity(Spin M1, Spin M2, std::string& nameDec, EvtData* theData);
  
  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);  
  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar);

  virtual void print(std::ostream& os) const;


protected:
  //  std::shared_ptr<resReaction> _resReactionPtr;
  std::shared_ptr<AbsXdecAmp> _XdecAmp;
  std::shared_ptr<PiPiScatteringXdecAmps> _PiPiScatteringXdecAmp;

  virtual void ThreadfuncData(unsigned int minEvent, unsigned int maxEvent,
                              std::shared_ptr<AbsPawianParameters> fitPar, LHData& theLHData);
  virtual void ThreadfuncMc(unsigned int minEvent, unsigned int maxEvent,
                            std::shared_ptr<AbsPawianParameters> fitPar, LHData& theLHData);
private:

  void initialize();
};
