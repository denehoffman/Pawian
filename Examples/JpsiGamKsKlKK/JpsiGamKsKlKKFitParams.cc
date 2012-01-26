#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKFitParams.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKStates.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

JpsiGamKsKlKKFitParams::JpsiGamKsKlKKFitParams() 
{
  JpsiGamKsKlKKStates theStates;
  theStates.print(std::cout);
  
  _jpclsMap[paramEnumJpsiGamKsKlKK::PsiToEtacGamma]=theStates.PsiToEtaGammaStates();
  _jpclsMap[paramEnumJpsiGamKsKlKK::PsiToEta2225Gamma]=theStates.PsiToEtaGammaStates();
  _jpclsMap[paramEnumJpsiGamKsKlKK::PsiToF02020Gamma]=theStates.PsiToF0GammaStates();
  _jpclsMap[paramEnumJpsiGamKsKlKK::PsiToF22010Gamma]=theStates.PsiToF2GammaStates();
  _jpclsMap[paramEnumJpsiGamKsKlKK::PsiToF22300Gamma]=theStates.PsiToF2GammaStates();
  _jpclsMap[paramEnumJpsiGamKsKlKK::PsiToF22340Gamma]=theStates.PsiToF2GammaStates(); 
  _jpclsMap[paramEnumJpsiGamKsKlKK::PsiToEta21870Gamma]=theStates.PsiToEta2GammaStates();
  _jpclsMap[paramEnumJpsiGamKsKlKK::PsiToF1Gamma]=theStates.PsiToF1GammaStates();
  _jpclsMap[paramEnumJpsiGamKsKlKK::F02020ToPhiPhi]=theStates.F0ToPhiPhiStates();
  _jpclsMap[paramEnumJpsiGamKsKlKK::F22300ToPhiPhi]=theStates.F2ToPhiPhiStates();
  _jpclsMap[paramEnumJpsiGamKsKlKK::Eta21870ToPhiPhi]=theStates.Eta2ToPhiPhiStates();
  _jpclsMap[paramEnumJpsiGamKsKlKK::F1ToPhiPhi]=theStates.F1ToPhiPhiStates();
}

JpsiGamKsKlKKFitParams::JpsiGamKsKlKKFitParams(fitParams& theStartparams, fitParams& theErrorparams) :
  FitParamsBase(theStartparams, theErrorparams)
{
  JpsiGamKsKlKKStates theStates;
  theStates.print(std::cout);
  
  _jpclsMap[paramEnumJpsiGamKsKlKK::PsiToEtacGamma]=theStates.PsiToEtaGammaStates();
  _jpclsMap[paramEnumJpsiGamKsKlKK::PsiToEta2225Gamma]=theStates.PsiToEtaGammaStates();
  _jpclsMap[paramEnumJpsiGamKsKlKK::PsiToF02020Gamma]=theStates.PsiToF0GammaStates();
  _jpclsMap[paramEnumJpsiGamKsKlKK::PsiToF22010Gamma]=theStates.PsiToF2GammaStates();
  _jpclsMap[paramEnumJpsiGamKsKlKK::PsiToF22300Gamma]=theStates.PsiToF2GammaStates();
  _jpclsMap[paramEnumJpsiGamKsKlKK::PsiToF22340Gamma]=theStates.PsiToF2GammaStates(); 
  _jpclsMap[paramEnumJpsiGamKsKlKK::PsiToEta21870Gamma]=theStates.PsiToEta2GammaStates(); 
  _jpclsMap[paramEnumJpsiGamKsKlKK::PsiToF1Gamma]=theStates.PsiToF1GammaStates();
  _jpclsMap[paramEnumJpsiGamKsKlKK::F02020ToPhiPhi]=theStates.F0ToPhiPhiStates();
  _jpclsMap[paramEnumJpsiGamKsKlKK::F22300ToPhiPhi]=theStates.F2ToPhiPhiStates();
  _jpclsMap[paramEnumJpsiGamKsKlKK::Eta21870ToPhiPhi]=theStates.Eta2ToPhiPhiStates();
  _jpclsMap[paramEnumJpsiGamKsKlKK::F1ToPhiPhi]=theStates.F1ToPhiPhiStates();
  
}







JpsiGamKsKlKKFitParams::~JpsiGamKsKlKKFitParams()
{
}



