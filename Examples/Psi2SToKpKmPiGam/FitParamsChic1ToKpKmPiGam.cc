#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2SToKpKmPiGam/FitParamsChic1ToKpKmPiGam.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamStates.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

FitParamsChic1ToKpKmPiGam::FitParamsChic1ToKpKmPiGam() :
  FitParamsBase() 
{
  Psi2SToKpKmPiGamStates theStates;
  theStates.print(std::cout);

  _jpclsMap[paramEnumChic1ToKpKmPiGam::ChiGam]=theStates.PsiToChiGamStates();
  _jpclsMap[paramEnumChic1ToKpKmPiGam::K890K]=theStates.ChiToKst1KStates();
  _jpclsMap[paramEnumChic1ToKpKmPiGam::K_0_1400K]=theStates.ChiToKst0KStates();
  _jpclsMap[paramEnumChic1ToKpKmPiGam::K_1_1400K]=theStates.ChiToKst1KStates();
  _jpclsMap[paramEnumChic1ToKpKmPiGam::K_2_1400K]=theStates.ChiToKst2KStates();
  _jpclsMap[paramEnumChic1ToKpKmPiGam::K_1_1680K]=theStates.ChiToKst1KStates();
  _jpclsMap[paramEnumChic1ToKpKmPiGam::a980Pi]=theStates.ChiToa0PiStates();
  _jpclsMap[paramEnumChic1ToKpKmPiGam::KappaK]=theStates.ChiToKst0KStates();
  _jpclsMap[paramEnumChic1ToKpKmPiGam::KKPi]=theStates.ChiToa0PiStates();
}

FitParamsChic1ToKpKmPiGam::FitParamsChic1ToKpKmPiGam(fitParams& theStartparams, fitParams& theErrorparams) :
  FitParamsBase(theStartparams, theErrorparams) 
{
  Psi2SToKpKmPiGamStates theStates;
  theStates.print(std::cout);

  _jpclsMap[paramEnumChic1ToKpKmPiGam::ChiGam]=theStates.PsiToChiGamStates();
  _jpclsMap[paramEnumChic1ToKpKmPiGam::K890K]=theStates.ChiToKst1KStates();
  _jpclsMap[paramEnumChic1ToKpKmPiGam::K_0_1400K]=theStates.ChiToKst0KStates();
  _jpclsMap[paramEnumChic1ToKpKmPiGam::K_1_1400K]=theStates.ChiToKst1KStates();
  _jpclsMap[paramEnumChic1ToKpKmPiGam::K_2_1400K]=theStates.ChiToKst2KStates();
  _jpclsMap[paramEnumChic1ToKpKmPiGam::K_1_1680K]=theStates.ChiToKst1KStates();
  _jpclsMap[paramEnumChic1ToKpKmPiGam::a980Pi]=theStates.ChiToa0PiStates();
  _jpclsMap[paramEnumChic1ToKpKmPiGam::KappaK]=theStates.ChiToKst0KStates();
  _jpclsMap[paramEnumChic1ToKpKmPiGam::KKPi]=theStates.ChiToa0PiStates();

}

FitParamsChic1ToKpKmPiGam::~FitParamsChic1ToKpKmPiGam()
{
}

