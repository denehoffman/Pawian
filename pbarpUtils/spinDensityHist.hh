// pbarpEventList class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#pragma once

#include <map>
#include <vector>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh" 

#include <boost/shared_ptr.hpp> 

class AbsLh;
class fitParams;
class EvtData;

class TH1F;
class TFile;

typedef std::map<int, std::map<int, TH1F*> > spinDensityHistMatrix;




class spinDensityHist{

  public:
    spinDensityHist(boost::shared_ptr<AbsLh> theLh, fitParams& theFitParams);
    ~spinDensityHist();
  private:
    short _nBins;
    int _maxEvents;
    TFile* _spinDensityRootFile;
    std::vector<EvtData*> _dataList;
    boost::shared_ptr<AbsLh> _theLh;

    void calcSpinDensityMatrix(std::string& particleName);
    void calcSpinDensityMatrixElement(std::string& particleName, Spin M1, Spin M2);
    void fillHistogram(std::string& particleName, TH1F* theHisto, EvtData* theData, double spinDensityValue);
};
