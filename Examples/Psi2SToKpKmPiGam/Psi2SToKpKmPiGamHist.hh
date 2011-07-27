#ifndef _Psi2SToKpKmPiGamHist_H
#define _Psi2SToKpKmPiGamHist_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cassert>

#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamData.hh"
#include "Examples/Psi2SToKpKmPiGam/FitParamsKpKmPiGam.hh"

using Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData;

class Psi2SToKpKmPiGamEventList;
class AbsPsi2SToKpKmPiGamLh;
class TFile;
class TH2F;
class TH1F;
class TNtuple;

class Psi2SToKpKmPiGamHist {

public:

  // create/copy/destroy:

  ///Constructor 
  Psi2SToKpKmPiGamHist(boost::shared_ptr<const Psi2SToKpKmPiGamEventList>);
  Psi2SToKpKmPiGamHist(boost::shared_ptr<AbsPsi2SToKpKmPiGamLh>, paramKpKmPiGam&);

  /** Destructor */
  virtual ~Psi2SToKpKmPiGamHist();

  // Getters:
 
protected:


private:
  TFile* _theTFile;
  TH2F* _dalitzDataHist;
  TH2F* _dalitzMcHist;
  TH2F* _dalitzFittedHist;
  TH1F* _cosPsiDataHist;
  TH1F* _cosPsiMcHist;
  TH1F* _cosPsiFittedHist;
  TH1F* _cosChic1Toa0PiDataHist;
  TH1F* _cosChic1Toa0PiMcHist;
  TH1F* _cosChic1Toa0PiFittedHist;
  TH1F* _cosChic1ToK890KDataHist;
  TH1F* _cosChic1ToK890KMcHist;
  TH1F* _cosChic1ToK890KFittedHist;
  TH1F* _cosChic1ToK1400KDataHist;
  TH1F* _cosChic1ToK1400KMcHist;
  TH1F* _cosChic1ToK1400KFittedHist;

  TH1F* _phiChic1ToK1400KDataHist;
  TH1F* _phiChic1ToK1400KMcHist;
  TH1F* _phiChic1ToK1400KFittedHist;

  TH1F* _cosChic1ToKKPiDataHist;
  TH1F* _cosChic1ToKKPiMcHist;
  TH1F* _cosChic1ToKKPiFittedHist;

  TH1F* _phiChic1ToKKPiDataHist;
  TH1F* _phiChic1ToKKPiMcHist;
  TH1F* _phiChic1ToKKPiFittedHist;

  TH1F* _cosK890DataHist;
  TH1F* _cosK890McHist;
  TH1F* _cosK890FittedHist;
  TH1F* _cosK1400DataHist;
  TH1F* _cosK1400McHist;
  TH1F* _cosK1400FittedHist;
  TH1F* _cosa980DataHist;
  TH1F* _cosa980McHist;
  TH1F* _cosa980FittedHist;
  TH1F* _invKpKmDataHist;
  TH1F* _invKpKmMcHist;
  TH1F* _invKpKmFittedHist;
  TH1F* _invKPiDataHist;
  TH1F* _invKPiMcHist;
  TH1F* _invKPiFittedHist; 
  TNtuple* _dataTuple;
  TNtuple* _mcTuple;



  void initRootStuff();
  void plotDalitz(TH2F* theHisto, const Psi2SToKpKmPiGamEvtData* theData, double weight);
  void plotCosPsi(TH1F* theHisto, const Psi2SToKpKmPiGamEvtData* theData, double weight);
  void plotChic1Toa0Pi(TH1F* theHisto, const Psi2SToKpKmPiGamEvtData* theData, double weight);
  void plotChic1ToK890K(TH1F* theHisto, const Psi2SToKpKmPiGamEvtData* theData, double weight);
  void plotChic1ToK1400K(TH1F* costHist, TH1F* phiHist,const Psi2SToKpKmPiGamEvtData* theData, double weight);

  void plotAnglesChic1ToKKPi(TH1F* costHist, TH1F* phiHist, const Psi2SToKpKmPiGamEvtData* theData, double weight);

  void plotCosK1400(TH1F* costHist, const Psi2SToKpKmPiGamEvtData* theData, double weight);
  void plotCosK890(TH1F* theHisto, const Psi2SToKpKmPiGamEvtData* theData, double weight);
  void plotCosa980(TH1F* theHisto, const Psi2SToKpKmPiGamEvtData* theData, double weight);
  void plotMKpKm(TH1F* theHisto, const Psi2SToKpKmPiGamEvtData* theData, double weight);
  void plotMKPi(TH1F* theHisto, const Psi2SToKpKmPiGamEvtData* theData, double weight);
  void fillTuple( TNtuple* theTuple, const Psi2SToKpKmPiGamEvtData* theData, double weight);

};

#endif
