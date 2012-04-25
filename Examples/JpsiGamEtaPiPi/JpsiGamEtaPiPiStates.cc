#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiStates.hh"
#include "ErrLogger/ErrLogger.hh"

JpsiGamEtaPiPiStates::JpsiGamEtaPiPiStates() : 
  JpsiToXGamStates(),
  _etaJPC(new jpcRes(0, -1, 1)),
  _f0JPC(new jpcRes(0, 1, 1)),
  _f1JPC(new jpcRes(1, 1, 1)),
  _f2JPC(new jpcRes(2, 1, 1)),
  _a0JPC(new jpcRes(0, 1, 1)),
  _piJPC(new jpcRes(0, -1, 1))
{

  //X decays to phi phi
  fillJPCLS(_etaJPC, _a0JPC, _piJPC, _JPCLS_EtaToa0Pi);

}

JpsiGamEtaPiPiStates::~JpsiGamEtaPiPiStates()
{
}



void JpsiGamEtaPiPiStates::print(std::ostream& os) const
{
  JpsiToXGamStates::print(os);
  os << "*** eta -> a0 pi:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_EtaToa0Pi );
  
}


