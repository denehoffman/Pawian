#ifndef _etacToa1320pi0fit_H
#define _etacToa1320pi0fit_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cassert>
#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/qft++/EtacToa1320pi0Data.hh"
#include "Minuit2/MnUserParameters.h"
//#include "Minuit2/FCNBase.h"

// using namespace std;
using namespace ROOT::Minuit2;

class TFile;
class TH2F;
class TH1F;

class EtacToa1320pi0fit {

public:

  // create/copy/destroy:

  ///Constructor 
  EtacToa1320pi0fit(int);



  /** Destructor */
  virtual ~EtacToa1320pi0fit();

  // Getters:

  double calcLogLh(const fitParamVal& theParamVal);
  bool initFitParameters(MnUserParameters& minuitParams);
  bool setFitParamVal(fitParamVal& fitParamVal, const std::vector<double>& par);
  bool fillFitHists(const fitParamVal& fitParamVal);
 
protected:


private:
  std::string _dataPath;
  std::string _mcPath;
  int _kindOfData;
  int _nOfData;
  int _nOfMc;
  bool _fitHistsFilled;


  TFile* _theTFile;
  TH2F* _dalitzDataHist;
  TH2F* _dalitzMcHist;
  TH2F* _dalitzFittedHist;
  TH1F* _invpietaDataHist;
  TH1F* _invpietaMcHist;
  TH1F* _invpietaFittedHist;
  TH1F* _invpipiDataHist;
  TH1F* _invpipiMcHist;
  TH1F* _invpipiFittedHist;

 std::vector<evt4Vec> _data4Vecs;
 std::vector<evt4Vec> _mc4Vecs;

  void initRootStuff();
  void read4Vecs(std::string& path,  int nEvts, std::vector<evt4Vec>& the4Vecs);
  
  void get4Vecs (std::ifstream& inStream, Vector4<double>& the4Vec);

  void plotDalitz (TH2F* theHisto, evt4Vec& theEvt4Vecs, double weight);
  void plotInvPiEta (TH1F* theHisto, evt4Vec& theEvt4Vecs, double weight);
  void plotInvPiPi (TH1F* theHisto, evt4Vec& theEvt4Vecs, double weight);

  double calcIntensityCache(evt4Vec& theEvtVec, const fitParamVal& theParamVal);

  Tensor<complex<double> > calcSpin2Amp(Vector4<double>& pi_cm, Vector4<double>& eta_cm, Vector4<double>& pi_recoil_cm);
  Tensor<complex<double> > calcSpin1Amp(Vector4<double>& pi_cm, Vector4<double>& eta_cm, Vector4<double>& pi_recoil_cm);

};

#endif /* _EtacToa1320pi0fit_H */
