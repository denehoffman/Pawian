#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiFitParams.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiStates.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

JpsiGamEtaPiPiFitParams::JpsiGamEtaPiPiFitParams()
{
  JpsiGamEtaPiPiStates theStates;
  theStates.print(std::cout);
  
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToEta1405Gamma]=theStates.PsiToEtaGammaStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToEta1295Gamma]=theStates.PsiToEtaGammaStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToF11285Gamma]=theStates.PsiToF1GammaStates();
  /*
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToF22010Gamma]=theStates.PsiToF2GammaStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToF22300Gamma]=theStates.PsiToF2GammaStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToF22340Gamma]=theStates.PsiToF2GammaStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToF1Gamma]=theStates.PsiToF1GammaStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToEta21870Gamma]=theStates.PsiToEta2GammaStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::F02020ToPhiPhi]=theStates.F0ToPhiPhiStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::F22300ToPhiPhi]=theStates.F2ToPhiPhiStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::Eta21870ToPhiPhi]=theStates.Eta2ToPhiPhiStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::F1ToPhiPhi]=theStates.F1ToPhiPhiStates();
  */
}

JpsiGamEtaPiPiFitParams::JpsiGamEtaPiPiFitParams(fitParams& theStartparams, fitParams& theErrorparams) :
  FitParamsBase(theStartparams, theErrorparams)
{
  JpsiGamEtaPiPiStates theStates;
  theStates.print(std::cout);
  
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToEta1405Gamma]=theStates.PsiToEtaGammaStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToEta1295Gamma]=theStates.PsiToEtaGammaStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToF11285Gamma]=theStates.PsiToF1GammaStates();
  /*
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToF22010Gamma]=theStates.PsiToF2GammaStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToF22300Gamma]=theStates.PsiToF2GammaStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToF22340Gamma]=theStates.PsiToF2GammaStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToF1Gamma]=theStates.PsiToF1GammaStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToEta21870Gamma]=theStates.PsiToEta2GammaStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::F02020ToPhiPhi]=theStates.F0ToPhiPhiStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::F22300ToPhiPhi]=theStates.F2ToPhiPhiStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::Eta21870ToPhiPhi]=theStates.Eta2ToPhiPhiStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::F1ToPhiPhi]=theStates.F1ToPhiPhiStates();
  */
}

JpsiGamEtaPiPiFitParams::~JpsiGamEtaPiPiFitParams()
{
}



