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

// gamgamStates class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <memory>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/AbsStates.hh"
#include "PwaUtils/DataUtils.hh"

class Particle;

class gamgamStates : public AbsStates {

public:

  gamgamStates();
  gamgamStates(int lmax);
  virtual ~gamgamStates();

  virtual std::vector< std::shared_ptr<const IGJPC> > igjpcStates(){
    return _allIGjpcRes; 
  }

  virtual std::vector< std::shared_ptr<const jpcRes> > jpcStates() const {
    return _alljpcRes;
  }

  virtual std::vector< std::shared_ptr<const JPCLS> > jpcLSStates() const {
    return _JPCLS_States;
  }

  virtual void print(std::ostream& os) const;

 protected:
  int _jmax;
  jpcRes _gammapJPC;
  jpcRes _gammaJPC;
  Particle* _gammaParticle;

  virtual bool calcStates();

  std::vector< std::shared_ptr<const IGJPC> > _allIGjpcRes;

  std::vector< std::shared_ptr<const jpcRes> > _alljpcRes;

  std::vector< std::shared_ptr<const JPCLS> > _JPCLS_States; //l;momentum between gamma and gamma
  std::vector< std::shared_ptr<const JPClamlam> > _JPClamlam_States;
                                                   
};
