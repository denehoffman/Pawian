#ifndef _Psi2STo2K2PiGamHist_H
#define _Psi2STo2K2PiGamHist_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cassert>

#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamData.hh"

using Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData;

class Psi2STo2K2PiGamEvtList;
class AbsPsi2STo2K2PiGamLh;
class TFile;
class TH2F;
class TH1F;
class TNtuple;

class Psi2STo2K2PiGamHist {

public:

  // create/copy/destroy:

  ///Constructor 
  Psi2STo2K2PiGamHist(boost::shared_ptr<const Psi2STo2K2PiGamEvtList>);
  Psi2STo2K2PiGamHist(boost::shared_ptr<AbsPsi2STo2K2PiGamLh>, Psi2STo2K2PiGamData::fitParamVal&);

  /** Destructor */
  virtual ~Psi2STo2K2PiGamHist();

  // Getters:
 
protected:


private:
  TFile* _theTFile;
  TH1F* _cosPsiDataHist;
  TH1F* _cosPsiMcHist;
  TH1F* _cosPsiFittedHist;
  TH1F* _cosK892DataHist;
  TH1F* _cosK892McHist;
  TH1F* _cosK892FittedHist;
  TH1F* _cosK1430DataHist;
  TH1F* _cosK1430McHist;
  TH1F* _cosK1430FittedHist;
  TH1F* _cosK1430ViaK892DataHist;
  TH1F* _cosK1430ViaK892McHist;
  TH1F* _cosK1430ViaK892FittedHist;
  TH1F* _invKKDataHist;
  TH1F* _invKKMcHist;
  TH1F* _invKKFittedHist;
  TH1F* _invKPiDataHist;
  TH1F* _invKPiMcHist;
  TH1F* _invKPiFittedHist;
  TH1F* _invPiPiDataHist;
  TH1F* _invPiPiMcHist;
  TH1F* _invPiPiFittedHist;

  TH1F* _invKPiPiViaK892DataHist;
  TH1F* _invKPiPiViaK892McHist;
  TH1F* _invKPiPiViaK892FittedHist;

  TH2F* _KPivsKPiDataHist2d;
  TH2F* _KPivsKPiMcHist2d;
  TH2F* _KPivsKPiFittedHist2d;
  TH2F* _KKvsPiPiDataHist2d;
  TH2F* _KKvsPiPiMcHist2d;
  TH2F* _KKvsPiPiFittedHist2d;
  TH2F* _KPiPivsPiPiDataHist2d;
  TH2F* _KPiPivsPiPiMcHist2d;
  TH2F* _KPiPivsPiPiFittedHist2d;
  TH2F* _KKPivsKPiDataHist2d;
  TH2F* _KKPivsKPiMcHist2d;
  TH2F* _KKPivsKPiFittedHist2d;

  TNtuple* _dataTuple;
  TNtuple* _mcTuple;

  void initRootStuff();
  void writeNTuple(TNtuple* theTuple, const Psi2STo2K2PiGamEvtData* theData, double weight);
  void plotCosPsi(TH1F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight);
  void plotCosKst(TH1F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight, double mass, double width);
  void plotCosKstViaK892(TH1F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight, double mass, double width);
  void plotInvKK(TH1F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight);
  void plotInvKPi(TH1F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight);
  void plotInvPiPi(TH1F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight);
  void plotInvKPiPi(TH1F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight, double mass, double width);
  void plotKPivsKPi(TH2F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight);
  void plotKKvsPiPi(TH2F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight);
  void plotKPiPivsPiPi(TH2F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight);
  void plotKKPivsKPi(TH2F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight);
  void fillTuple( TNtuple* theTuple, const Psi2STo2K2PiGamEvtData* theData, double weight);

};

#endif
