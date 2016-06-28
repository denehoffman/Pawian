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

#ifndef _BwShape_H
#define _BwShape_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cassert>

#include <memory>

#include "TROOT.h"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

class TFile;
class TH1F;
class TH2F;

class BwShape {

public:

  // create/copy/destroy:

  ///Constructor 
  BwShape(double MassRes, double MassWidth, double MassDec1, double MassDec2, unsigned int Lmax, double deltaMass);


  /** Destructor */
  virtual ~BwShape();

  // Getters:
 
protected:


private:
  TFile* _theTFile;
  // std::map <unsigned int, TH1F* > _histMap;
  std::map <unsigned int, TH1F* > _histMapNew;
  // std::map <unsigned int, TH2F* > _argandHistMap;
  std::map <unsigned int, TH2F* > _argandHistMapNew;
  // std::map <unsigned int, TH2F* > _phaseHistMap;
};

#endif
