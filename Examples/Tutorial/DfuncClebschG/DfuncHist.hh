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

#ifndef _DfuncHist_H
#define _DfuncHist_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cassert>

#include <memory>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

class TH1F;
class TFile;

class DfuncHist {

public:

  // create/copy/destroy:

  ///Constructor 
  DfuncHist(Spin J, Spin lam1, Spin lam2);


  /** Destructor */
  virtual ~DfuncHist();

  // Getters:
  void print(std::ostream& os) const; 

protected:


private:
  Spin _J;
  Spin _lam1;
  Spin _lam2;
  TFile* _theTFile;
  TH1F* _dcostHist;
};

#endif
