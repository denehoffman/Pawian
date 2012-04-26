#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiFitParams.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiStates.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

JpsiGamEtaPiPiFitParams::JpsiGamEtaPiPiFitParams()
{
  initialize();
}

JpsiGamEtaPiPiFitParams::JpsiGamEtaPiPiFitParams(fitParams& theStartparams, fitParams& theErrorparams) :
  FitParamsBase(theStartparams, theErrorparams)
{
  initialize();
}

JpsiGamEtaPiPiFitParams::~JpsiGamEtaPiPiFitParams()
{
}

void JpsiGamEtaPiPiFitParams::initialize(){
  JpsiGamEtaPiPiStates theStates;
  theStates.print(std::cout);
  
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToEtaGamma]=theStates.PsiToEtaGammaStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToEta2Gamma]=theStates.PsiToEta2GammaStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToF1Gamma]=theStates.PsiToF0GammaStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::EtaToPiPiEta]=theStates.EtaToa0PiStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::EtaToA980Pi]=theStates.EtaToa0PiStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::EtaToA2_1320Pi]=theStates.EtaToa2PiStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::Eta2ToPiPiEta]=theStates.Eta2Toa0PiStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::Eta2ToA980Pi]=theStates.Eta2Toa0PiStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::Eta2ToA2_1320Pi]=theStates.Eta2Toa2PiStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::F1ToPiPiEta]=theStates.F1Toa0PiStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::F1ToA980Pi]=theStates.F1Toa0PiStates();
}

