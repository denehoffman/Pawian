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
  
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToEtaGamma]=theStates.PsiToEtaGammaStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToF1Gamma]=theStates.PsiToF0GammaStates();

}

JpsiGamEtaPiPiFitParams::JpsiGamEtaPiPiFitParams(fitParams& theStartparams, fitParams& theErrorparams) :
  FitParamsBase(theStartparams, theErrorparams)
{
  JpsiGamEtaPiPiStates theStates;
  theStates.print(std::cout);

  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToEtaGamma]=theStates.PsiToEtaGammaStates();
  _jpclsMap[paramEnumJpsiGamEtaPiPi::PsiToF1Gamma]=theStates.PsiToF1GammaStates();  
}

JpsiGamEtaPiPiFitParams::~JpsiGamEtaPiPiFitParams()
{
}



