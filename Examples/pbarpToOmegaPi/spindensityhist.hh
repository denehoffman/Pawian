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
#ifndef SPINDENSITYHIST_H
#define SPINDENSITYHIST_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cassert>

#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiData.hh"

using OmegaPiData::OmPiEvtData;

class OmegaPiEventList;
class OmegaPiLh;
class TFile;
class TH2F;
class TH1F;


/**
Spin Density Histogrammer. It is used to create histograms from calculated spin density.

	@author Denis Pavlina <denis@pc15>
*/
class SpinDensityHist
{
  public:
    ~SpinDensityHist();
    
    SpinDensityHist(const std::string &thePathToRootFile, const std::vector<OmegaPiData::OmPiEvtData*> &theEventData,
                    OmegaPiData::fitParamVal &theParamVal);
    void createHistogram(int M, int M_, bool bIncoherent);
    void createHistograms();

  protected:
    unsigned int _numOfEvts;
    const std::vector<OmegaPiData::OmPiEvtData*> m_EventData;
    OmegaPiData::fitParamVal *m_pfitParamVal;
    TFile* _theTFile;
    TH1F* _cosOmegaHeliMcHist;
    std::string m_PathToRootFile;
    void createSpinDensityHist(TH1F* theHisto, int M, int M_, bool bIncoherent,bool bReal=true);

    void AddData(TH1F* theHisto, const OmPiEvtData& theEvtData, double dSpinDensity);
};

#endif
