//************************************************************************//
//									  //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// pbarpState class definition file. -*- C++ -*-
// Copyright 2020 Bertram Kopf

#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <memory>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/DataUtils.hh"


class pbarpState {

public:

  pbarpState(std::shared_ptr<const jpcRes> theJPC);
  virtual ~pbarpState();

  std::shared_ptr<const jpcRes> jpcState() const {return _jpc;} 
  bool isSingletState() const {return _isSinglet;}
  bool isTripletState() const {return _isTriplet;}
  bool isTriplet0State() const {return _isTriplet0;}
  bool isTriplet1State() const {return _isTriplet1;}
 
  virtual void print(std::ostream& os) const;

 protected:
  std::shared_ptr<const jpcRes> _jpc;
  jpcRes _pbarJPC;
  jpcRes _pJPC;
  vector<LS> _pbarpLSs;
  bool _ispbarpState;
  bool _isSinglet;
  bool _isTriplet;
  bool _isTriplet0;
  bool _isTriplet1;
  
  void extractLSs();
};
