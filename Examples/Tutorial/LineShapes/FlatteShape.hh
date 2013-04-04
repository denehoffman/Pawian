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

#ifndef _FlatteShape_H
#define _FlatteShape_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cassert>

#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

class TFile;
class TH1F;
class TH2F;

class FlatteShape {

public:

  // create/copy/destroy:

  ///Constructor 
  FlatteShape(std::string ptype, double g1, double g2);


  /** Destructor */
  virtual ~FlatteShape();

  // Getters:
 
protected:


private:
  TFile* _theTFile;
  TH1F* _histShapeLow;
  TH1F* _histShapeHigh;
  TH1F* _histShapeLowFlatteNew;
  TH1F* _histShapeHighFlatteNew;
  TH1F* _histShapeLowTest;
  TH1F* _histShapeHighTest;
  TH1F* _histShapeLowKmatr;
  TH1F* _histShapeHighKmatr;
  TH2F* _argandKmatrFlatteLowHist;
  TH2F* _argandKmatrLowHist;
  TH2F* _argandKmatrHighHist;

  TH1F* _histShapeLowKmatrRel;
  TH1F* _histShapeHighKmatrRel;
  TH2F* _argandKmatrLowRelHist;
  TH2F* _argandKmatrHighRelHist;


  std::string _pType;

};

#endif
