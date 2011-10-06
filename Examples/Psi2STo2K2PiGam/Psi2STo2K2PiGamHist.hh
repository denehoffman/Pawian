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
#include "Examples/Psi2STo2K2PiGam/FitParams2K2PiGam.hh"

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
  Psi2STo2K2PiGamHist(boost::shared_ptr<const Psi2STo2K2PiGamEvtList>, std::string fileName="Psi2STo2K2PiGamQa.root");
  Psi2STo2K2PiGamHist(boost::shared_ptr<AbsPsi2STo2K2PiGamLh>, param2K2PiGam&, std::string fileName="Psi2STo2K2PiGamPWA.root");
  Psi2STo2K2PiGamHist(const param2K2PiGam&, std::string fileName="Psi2STo2K2PiGamGen.root");

  /** Destructor */
  virtual ~Psi2STo2K2PiGamHist();

  void fillHistos(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight);
  // Getters:
 
protected:


private:
  TFile* _theTFile;
  std::map <std::string, TH1F* > _hist1DMap;
  std::map <std::string, TH2F* > _hist2DMap;
  const std::string _dataKey;
  const std::string _fittedKey;
  const std::string _mcKey;
  std::vector<std::string> _keyVec;


  TNtuple* _dataTuple;
  TNtuple* _mcTuple;

  void initRootStuff(const std::string& fileName);
  void writeNTuple(TNtuple* theTuple, const Psi2STo2K2PiGamEvtData* theData, double weight);
  void plotCosPsi(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight);
  void plotCosKst(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight, double mass, double width);
  void plotCosKstViaK892(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight, double mass, double width);
  void plotInvKK(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight);
  void plotInvKPi(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight);
  void plotInvPiPi(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight);
  void plotInvKPiPi(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight, double mass, double width);
  void plotKPivsKPi(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight);
  void plotKKvsPiPi(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight);
  void plotKPiPivsPiPi(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight);
  void plotKKPivsKPi(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight);


  void fillTuple( TNtuple* theTuple, const Psi2STo2K2PiGamEvtData* theData, double weight);

  void initHistMap();

};

#endif
