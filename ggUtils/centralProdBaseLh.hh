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

// centralProdBaseLh class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf, Malte Albrecht

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
#include "ggUtils/GGChannelEnv.hh"

class AbsXdecAmp;
class ggReaction;
class LSDecAmps;

class centralProdBaseLh : public AbsLh {

public:
  centralProdBaseLh(ChannelID channelID);

  virtual ~centralProdBaseLh();

  virtual AbsLh* clone_() const{
    AbsLh* theClone=new centralProdBaseLh(_channelID);
    theClone->setDataVec(_evtDataVec);
    theClone->setMcVec(_evtMCVec);
    return theClone;
  }

  virtual double calcEvtIntensity( EvtData* theData, std::shared_ptr<AbsPawianParameters> fitPar);
  virtual complex<double> calcProdPartAmp(Spin lamX, Spin lamDec, std::string nameDec, EvtData* theData,
					  std::map <std::shared_ptr<const JPCLS>,
					  std::vector< std::shared_ptr<AbsXdecAmp> >,
					  pawian::Collection::SharedPtrLess > pbarpAmps);

  virtual complex<double> calcSpinDensity(Spin M1, Spin M2, std::string& nameDec, EvtData* theData);

  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void print(std::ostream& os) const;


protected:
  const std::shared_ptr<GGChannelEnv> _ggChannelEnv;
  bool _useProdDynamics;
  std::vector<Spin> _lamX;
  std::shared_ptr<AbsDynamics> _dyn;

  virtual void initialize();
private:
};
