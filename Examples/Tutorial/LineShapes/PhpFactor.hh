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
  
  TH2F* _CMDefaultComplRealHist2;
  TH2F* _CMDefaultComplImagHist2;
  TH2F* _CMPenningtonRealHist2;
  TH2F* _CMPenningtonImagHist2;
  TH2F* _CMReidRealHist2;
  TH2F* _CMReidImagHist2;
  TH2F* _CMDudekRealHist2;
  TH2F* _CMDudekImagHist2;

  TH1F* _CMAsnerRealHist;
  TH1F* _CMAsnerImagHist;
  TH1F* _CMPenningtonRealHist;
  TH1F* _CMPenningtonImagHist;
  TH1F* _CMReidRealHist;
  TH1F* _CMReidImagHist;
  TH1F* _CMDudekRealHist;
  TH1F* _CMDudekImagHist;
  TH1F* _CMDefaultRealHist;
  TH1F* _CMDefaultImagHist;
  TH1F* _CMDefaultComplRealHist;
  TH1F* _CMDefaultComplImagHist;

  TH1F* _CMMAsnerRealHist;
  TH1F* _CMMAsnerImagHist;
  TH1F* _CMMPenningtonRealHist;
  TH1F* _CMMPenningtonImagHist;
  TH1F* _CMMReidRealHist;
  TH1F* _CMMReidImagHist;
  TH1F* _CMMDudekRealHist;
  TH1F* _CMMDudekImagHist;
  TH1F* _CMMDefaultComplRealHist;
  TH1F* _CMMDefaultComplImagHist;
  TH1F* _CMMBBUnstableRhoPiRealHist;
  TH1F* _CMMBBUnstableRhoPiImagHist;

  const double _mass1;
  const double _mass2;
  const double _massMax;
};


