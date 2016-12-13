//************************************************************************//
//									  //
//  Copyright 2016 Bertram Kopf (bertram@ep1.rub.de)			  //
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

#include <memory>

#include "TROOT.h"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

class TFile;
class TH1F;
class TH2F;

class PhpFactor {

public:

  // create/copy/destroy:

  ///Constructor 
  PhpFactor(double mass1, double mass2, double massMax);


  /** Destructor */
  virtual ~PhpFactor();

  // Getters:
 
protected:


private:
  TFile* _theTFile;
  TH1F* _phpRealHist;
  TH1F* _phpImagHist;

  TH2F* _phpDefaultComplRealHist2;
  TH2F* _phpDefaultComplImagHist2;
  TH2F* _phpPenningtonRealHist2;
  TH2F* _phpPenningtonImagHist2;
  TH2F* _phpReidRealHist2;
  TH2F* _phpReidImagHist2;

  TH1F* _phpAsnerRealHist;
  TH1F* _phpAsnerImagHist;
  TH1F* _phpPenningtonRealHist;
  TH1F* _phpPenningtonImagHist;
  TH1F* _phpReidRealHist;
  TH1F* _phpReidImagHist;
  TH1F* _phpDefaultRealHist;
  TH1F* _phpDefaultImagHist;
  TH1F* _phpDefaultComplRealHist;
  TH1F* _phpDefaultComplImagHist;

  TH1F* _phpMAsnerRealHist;
  TH1F* _phpMAsnerImagHist;
  TH1F* _phpMPenningtonRealHist;
  TH1F* _phpMPenningtonImagHist;
  TH1F* _phpMReidRealHist;
  TH1F* _phpMReidImagHist;
  TH1F* _phpMDefaultComplRealHist;
  TH1F* _phpMDefaultComplImagHist;

  const double _mass1;
  const double _mass2;
  const double _massMax;
};


