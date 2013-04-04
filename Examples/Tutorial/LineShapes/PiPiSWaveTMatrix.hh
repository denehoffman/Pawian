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
//PiPiSWaveTMatrix class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

class TFile;
class TH1F;
class TH2F;


class PiPiSWaveTMatrix {

public:

  // create/copy/destroy:

  ///Constructor 
  PiPiSWaveTMatrix();


  /** Destructor */
  virtual ~PiPiSWaveTMatrix();

  // Getters:
 
protected:


private:
  TFile* _theTFile;
  TH1F* _invPiPiMassH1;
  TH1F* _invPiPiMassRelH1;
  TH1F* _absT00RelH1;
  TH1F* _absS00RelH1;
  TH1F* _sqrT00RelSigmaPoleH1;
  TH2F* _pipiPhaseSpaceFactorH2;
  TH2F* _pipipipiPhaseSpaceFactorH2;
  TH2F*  _argandH2;
  TH2F*  _argandRelH2;
  TH2F* _phaseShiftH2;
  TH2F* _phaseShiftRelH2;
};


