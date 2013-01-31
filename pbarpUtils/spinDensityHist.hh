// pbarpEventList class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#pragma once

#include <map>
#include <vector>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh" 
#include <boost/shared_ptr.hpp> 
#include "Minuit2/MnUserParameters.h"
#include "PwaUtils/FitParamsBase.hh"

class AbsLh;
class fitParams;
class EvtData;
class PwaCovMatrix;

class TH1F;
class TFile;

typedef std::map<int, std::map<int, TH1F*> > spinDensityHistMatrix;




class spinDensityHist{

  public:
   spinDensityHist(boost::shared_ptr<AbsLh> theLh, fitParams& theFitParams, PwaCovMatrix& thePwaCovMatrix);
   ~spinDensityHist();

  private:
   short _nBins;
   int _maxEvents;
   TFile* _spinDensityRootFile;
   FitParamsBase theFitParamsBaseClass;
   std::vector<EvtData*> _dataList;
   boost::shared_ptr<AbsLh> _theLh;
   fitParams* _theFitParamsOriginal;
   PwaCovMatrix* _thePwaCovMatrix;
   ROOT::Minuit2::MnUserParameters _theMnUserParameters;

   void calcSpinDensityMatrix(std::string& particleName);
   void calcSpinDensityMatrixElement(std::string& particleName, Spin M1, Spin M2);
   complex<double> calcSpinDensityMatrixError(std::string& particleName, Spin M1, Spin M2, EvtData* evtData, complex<double> sdmValue);
   void fillHistogram(std::string& particleName, TH1F* theHisto, EvtData* theData, double spinDensityValue);
};
