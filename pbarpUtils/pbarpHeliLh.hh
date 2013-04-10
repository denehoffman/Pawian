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

// pbarpHeliLh class definition file. -*- C++ -*-
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

#include "pbarpUtils/pbarpBaseLh.hh"
#include "PwaUtils/DataUtils.hh"
#include "Minuit2/MnUserParameters.h"

class AbsXdecAmp;
class pbarpReaction;
class HeliDecAmps;

class pbarpHeliLh : public pbarpBaseLh {

public:
  // pbarpHeliLh(boost::shared_ptr<const EvtDataBaseList>);
  pbarpHeliLh();  
  virtual ~pbarpHeliLh();
  
  virtual AbsLh* clone_() const {
    AbsLh* theClone=new pbarpHeliLh();
    theClone->setDataVec(_evtDataVec);
    theClone->setMcVec(_evtMCVec);
    return theClone;
    //    return new  pbarpHeliLh(_evtListPtr);
  }
  virtual void print(std::ostream& os) const;
  
protected:
  
private:
  
  void initialize();
  
};
