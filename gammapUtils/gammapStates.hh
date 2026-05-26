//************************************************************************//
//                                                                        //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)                      //
//               Julian Pychy (julian@ep1.rub.de)                         //
//               - Ruhr-Universität Bochum                                //
//                                                                        //
//  This file is part of Pawian.                                          //
//                                                                        //
//  Pawian is free software: you can redistribute it and/or modify        //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or     //
//  (at your option) any later version.                                   //
//                                                                        //
//  Pawian is distributed in the hope that it will be useful,             //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of        //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         //
//  GNU General Public License for more details.                          //
//                                                                        //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.       //
//                                                                        //
//************************************************************************//

// gammapStates class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#pragma once

#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "PwaUtils/AbsStates.hh"
#include "PwaUtils/DataUtils.hh"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

class gammapStates : public AbsStates {

public:
  gammapStates();
  gammapStates(int lmax);
  virtual ~gammapStates();

  virtual std::vector<std::shared_ptr<const IGJPC>> igjpcStates() {
    return _allIGjpcRes;
  }

  virtual std::vector<std::shared_ptr<const jpcRes>> jpcStates() const {
    return _alljpcRes;
  }

  virtual std::vector<std::shared_ptr<const JPCLS>> jpclsStates() const {
    return _JPCLS_states;
  }

  virtual void print(std::ostream &os) const;

protected:
  int _lmax;
  jpcRes _gammapJPC;
  jpcRes _pJPC;
  jpcRes _gammaJPC;

  virtual bool calcStates();

  std::vector<std::shared_ptr<const IGJPC>> _allIGjpcRes;

  std::vector<std::shared_ptr<const jpcRes>> _alljpcRes;

  std::vector<std::shared_ptr<const JPCLS>>
      _JPCLS_states; // l;momentum between p and gamma
                     // s multipole gamma l coupling
  std::vector<std::shared_ptr<const JPCLS>> _elMpolJPCLS_States;
  std::vector<std::shared_ptr<const JPCLS>> _magMpolJPCLS_States;
};
