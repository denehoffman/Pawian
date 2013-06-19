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

#ifndef _TwoPoles_H
#define _TwoPoles_H

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

class TwoPoles {

public:

  // create/copy/destroy:

  ///Constructor 
  TwoPoles(double MassRes1, double Width1, double MassRes2, double Width2, double deltaMass);


  /** Destructor */
  virtual ~TwoPoles();

  // Getters:
 
protected:


private:
  TFile* _theTFile;
  TH1F* _massShapeHist;
  TH2F* _argandHist;
  TH2F* _phaseHist;
  TH1F* _massShapeKmatrHist;
  TH2F* _argandKmatrHist;
  TH2F* _phaseKmatrHist;
  TH1F* _massShapeKmatrHistRel;
  TH2F* _argandKmatrHistRel;
  TH2F* _phaseKmatrHistRel;
  double _massLow;
  double _widthLow;
  double _massHigh;
  double _widthHigh;
};

#endif
