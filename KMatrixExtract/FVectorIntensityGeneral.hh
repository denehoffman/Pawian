//************************************************************************//
//									  //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)			  //
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
//FVectorIntensityGeneral class definition file. -*- C++ -*-
// Copyright 2020 Bertram Kopf

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <memory>

#include "TROOT.h"
#include "KMatrixExtract/TMatrixGeneral.hh"

class TFile;
class TH1F;
class TH2F;
class KMatrixParser;
class AbsPhaseSpace;
class TMatrixRel;
class KMatrixRel;
class ParticleTable;
class KMatrixParser;
class AbsPawianParameters;
class pipiScatteringParser;
class FVectorIntensityDynamics;
class FVector;
class TGraph;
class PiPiScatteringChannelEnv;

class FVectorIntensityGeneral : public TMatrixGeneral{

public:

  // create/copy/destroy:

  ///Constructor 
  FVectorIntensityGeneral(pipiScatteringParser* theParser);


  /** Destructor */
  virtual ~FVectorIntensityGeneral();

  virtual void fillParams();
  virtual void process();
  virtual void initHistos();

  // Getters:
 
protected:


private:
  std::shared_ptr<FVectorIntensityDynamics> _fVectorIntensityDynamics;
  std::shared_ptr<FVector> _fVector; 
  std::string _pVecName;

  std::vector<TH1F*> _MagsH1;
  std::vector<TH1F*> _PhasesH1;
  std::vector<TH1F*> _IntensitiesH1;
  std::vector<TGraph*> _ArgandPlotsTGraph;  

  void init();
};


