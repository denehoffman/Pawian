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

#ifndef _JPCdecays_H
#define _JPCdecays_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cassert>

#include <memory>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

class jpcRes;

class JPCdecays {

public:

  // create/copy/destroy:

  ///Constructor 
  JPCdecays(std::shared_ptr<const jpcRes> motherJPC, std::shared_ptr<const jpcRes> daughter1JPC, std::shared_ptr<const jpcRes> daughter2JPC, bool sameDaughters=false);


  /** Destructor */
  virtual ~JPCdecays();

  // Getters:
  void print(std::ostream& os) const; 

protected:


private:
  std::shared_ptr<const jpcRes> _motherJPC;
  std::shared_ptr<const jpcRes> _daughter1JPC;
  std::shared_ptr<const jpcRes> _daughter2JPC;
  std::vector<LS> _allLSs;
  std::vector<LS> _allowedCParityLSs;
  std::vector<LS> _forbiddenCParityLSs;
  bool _sameDaughters;
};

#endif
