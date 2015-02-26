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

// spinDensityHist class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#pragma once

#include <map>
#include <vector>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh" 
#include <memory> 
//#include "Minuit2/MnUserParameters.h"
#include "PwaUtils/FitParamsBase.hh"

class AbsLh;
class fitParams;
class EvtData;
class PwaCovMatrix;
class AbsPawianParameters;

class TH1F;
class TFile;

typedef std::map<int, std::map<int, TH1F*> > spinDensityHistMatrix;




class spinDensityHist{

  public:
   spinDensityHist(std::shared_ptr<AbsLh> theLh, fitParams& theFitParams);
   ~spinDensityHist();
   void SetCovarianceMatrix(std::shared_ptr<PwaCovMatrix> thePwaCovMatrix);
   void Calculate();

   static const int MAX_EVENTS;

  private:
   bool _calcErrors;
   short _nBins;
   TFile* _spinDensityRootFile;
   FitParamsBase theFitParamsBaseClass;
   std::vector<EvtData*> _dataList;
   std::shared_ptr<AbsLh> _theLh;
   fitParams* _theFitParamsOriginal;
   std::shared_ptr<PwaCovMatrix> _thePwaCovMatrix;
  //   ROOT::Minuit2::MnUserParameters _theMnUserParameters;

  std::shared_ptr<AbsPawianParameters> _theParameters;

   void calcSpinDensityMatrix(std::string& particleName);
   void calcSpinDensityMatrixNorm(std::string& particleName, int J);
   void calcSpinDensityMatrixElement(std::string& particleName, Spin M1, Spin M2, int J);
   complex<double> calcSpinDensityMatrixError(std::string& particleName, Spin M1, Spin M2, EvtData* evtData, complex<double> sdmValue, int J);
   double ParticleCosTheta(std::string& particleName, EvtData* theData);
};
