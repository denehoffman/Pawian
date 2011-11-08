#ifndef _OmegaPiHistLS_H
#define _OmegaPiHistLS_H

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
#include "TH3F.h"
// #include "Examples/pbarpToOmegaPiLS/AbsOmegaPiEventListLS.hh"

using OmegaPiDataLS::OmPiEvtDataLS;

class AbsOmegaPiEventListLS;
class AbsOmegaPiLhLS;
class TFile;
class TH2F;
class TH1F;


class OmegaPiHistLS {

public:

  // create/copy/destroy:

  ///Constructor 
  OmegaPiHistLS(boost::shared_ptr<const AbsOmegaPiEventListLS>, const std::string &thePathToRootFile);
  OmegaPiHistLS(boost::shared_ptr<AbsOmegaPiLhLS>, OmegaPiDataLS::fitParamVal&, const std::string &thePathToRootFile);

  /** Destructor */
  virtual ~OmegaPiHistLS();

  // Getters:
 
protected:


private:
  TFile* _theTFile;
  TH1F* _cosOmegaHeliDataHist;
  TH1F* _cosOmegaHeliMcHist;
  TH1F* _cosOmegaHeliFittedHist;
  TH1F* _cosOmegaAccCorHist;
  TH1F* _cosOmegaFittedAccCorHist;

  TH1F* _cosPi0FromOmegaDataHeli;
  TH1F* _cosPi0FromOmegaMcHeli;
  TH1F* _cosPi0FromOmegaFittedHeli;
  TH1F* _cosPi0FromOmegaAccCorHeli;
  TH1F* _cosPi0FromOmegaFittedAccCorHeli;

  TH1F* _treimanYangDataHist;
  TH1F* _treimanYangMcHist;
  TH1F* _treimanYangFittedHist;
  TH1F* _treimanYangAccCorHist;
  TH1F* _treimanYangFittedAccCorHist;

  TH1F* _cosPi0FromOmegaDataHeli1; 

  TH1F* _lambdaOmegaTo3PiData;
  TH1F* _lambdaOmegaTo3PiMc;
  TH1F* _lambdaOmegaTo3PiFitted; 

  TH2F* _thetaPhiPi0FromOmegaDataHeli;
  TH2F* _thetaPhiPi0FromOmegaMcHeli;
  TH2F* _thetaPhiPi0FromOmegaFittedHeli;
  TH2F* _thetaPhiPi0FromOmegaAccCorHeli;
  TH2F* _thetaPhiPi0FromOmegaFittedAccCorHeli;

  TH3F* _prodDecThetaPhiPi0FromOmegaDataHeli;
  TH3F* _prodDecThetaPhiPi0FromOmegaMcHeli;
  TH3F* _prodDecThetaPhiPi0FromOmegaFittedHeli;

  unsigned _lmax;
  unsigned _pbarmom;

  void initRootStuff(const std::string &thePathToRootFile);
  void plotCosOmegaHeli(TH1F* theHisto, const OmPiEvtDataLS* theEvtData, double weight);
  void plotCosPi0FromOmegaHeli(TH1F* theHisto, const OmPiEvtDataLS* theEvtData, double weight);
  void plotCosPi0FromOmegaHeli1(TH1F* theHisto, const OmPiEvtDataLS* theEvtData, double weight);
  void plotTreimanYang(TH1F* theHisto, const OmPiEvtDataLS* theEvtData, double weight);
  void plotThetaPhiPi0FromOmegaHeli(TH2F* theHisto, const OmPiEvtDataLS* theEvtData, double weight);
  void plotProdDecThetaPhiPi0FromOmegaHeli(TH3F* theHisto, const OmPiEvtDataLS* theEvtData, double weight);
  void plotLambdaOmegaTo3Pi(TH1F* theHisto, const OmPiEvtDataLS* theEvtData, double weight);
  TH1F* doAccCor(TH1F* dataHist, TH1F* mcHist, const std::string& name, const std::string& title);
  TH2F* doAccCor(TH2F* dataHist, TH2F* mcHist, const std::string& name, const std::string& title);
};

#endif
