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


#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include <cassert>



#include "TROOT.h"
// #include <TSystem.h>
//KPiSWaveTMatrix class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

class TFile;
class TH1F;
class TH2F;


class KPiSWaveTMatrix {

public:

  // create/copy/destroy:

  ///Constructor 
  KPiSWaveTMatrix();


  /** Destructor */
  virtual ~KPiSWaveTMatrix();

  // Getters:
 
protected:


private:
  TFile* _theTFile;
  TH1F* _KPiAmpRealH1;
  TH1F* _KPiAmpImagH1;
  TH1F* _KPiAmpRealLASSH1;
  TH1F* _KPiAmpImagLASSH1;
  TH1F* _KPiArgandUnitH1;
  TH1F* _KPiPhaseH1;
};


