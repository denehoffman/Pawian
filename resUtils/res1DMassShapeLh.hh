//************************************************************************//
//									  //
//  Copyright 2021 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// res1DMassShapeLh class definition file. -*- C++ -*-
// Copyright 2021 Bertram Kopf

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <memory>
#include <boost/function.hpp>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

#include "resUtils/resBaseLh.hh"


class res1DMassShapeLh : public resBaseLh {

public:
  // resBaseLh(std::shared_ptr<const EvtDataBaseList>);
  res1DMassShapeLh(ChannelID channelID);

  virtual ~res1DMassShapeLh();

  virtual AbsLh* clone_() const{
    AbsLh* theClone=new res1DMassShapeLh(_channelID);
    theClone->setDataVec(_evtDataVec);
    theClone->setMcVec(_evtMCVec);
    return theClone;
  }

  virtual double calcEvtIntensity( EvtData* theData, std::shared_ptr<AbsPawianParameters> fitPar);
  
  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);  
  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual bool checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld) {return true;}
  
  virtual void print(std::ostream& os) const;


protected:
  std::vector< std::shared_ptr<AbsDynamics> > _dynVecs;

private:

};
