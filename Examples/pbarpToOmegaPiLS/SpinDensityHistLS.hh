//
// C++ Interface: spindensityhist
//
// Description: 
//
//
// Author: Denis Pavlina <denis@pc15>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SPINDENSITYHISTHH_H
#define SPINDENSITYHISTHH_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cassert>

#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/pbarpToOmegaPiLS/OmegaPiDataLS.hh"

using OmegaPiDataLS::OmPiEvtDataLS;

class AbsOmegaPiEventListLS;
class AbsOmegaPiLhLS;
class TFile;
class TH2F;
class TH1F;


/**
Spin Density Histogrammer. It is used to create histograms from calculated spin density.

	@author Denis Pavlina <denis@pc15>
*/
class SpinDensityHistLS
{
  public:
    ~SpinDensityHistLS();
    
  SpinDensityHistLS(const std::string &thePathToRootFile, boost::shared_ptr<AbsOmegaPiLhLS> absOmegaPiLh,
		  OmegaPiDataLS::fitParamVal &theParamVal);
    void createHistogram(int M, int M_);
    void createHistograms();

 protected:
  unsigned int _numOfEvts;
  boost::shared_ptr<AbsOmegaPiLhLS> _omegaPiLh;
  std::vector<OmegaPiDataLS::OmPiEvtDataLS*> m_EventData;
  OmegaPiDataLS::fitParamVal *m_pfitParamVal;
  TFile* _theTFile;
  TH1F* _cosOmegaHeliMcHist;
  std::string m_PathToRootFile;
  void createSpinDensityHist(TH1F* theHisto, int M, int M_,bool bReal=true);
  
  void AddData(TH1F* theHisto, const OmPiEvtDataLS& theEvtData, double dSpinDensity);

};

#endif
