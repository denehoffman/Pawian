#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiStates.hh"
#include "ErrLogger/ErrLogger.hh"

JpsiGamEtaPiPiStates::JpsiGamEtaPiPiStates() : 
  JpsiToXGamStates(),
  _etaJPC(new jpcRes(0, -1, 1)),
  _eta2JPC(new jpcRes(2, -1, 1)),
  _f0JPC(new jpcRes(0, 1, 1)),
  _f1JPC(new jpcRes(1, 1, 1)),
  _f2JPC(new jpcRes(2, 1, 1)),
  _a0JPC(new jpcRes(0, 1, 1)),
  _a2JPC(new jpcRes(2, 1, 1)),
  _piJPC(new jpcRes(0, -1, 1))
{

  //X decays to phi phi
  fillJPCLS(_etaJPC, _a0JPC, _piJPC, _JPCLS_EtaToa0Pi);
  fillJPCLS(_etaJPC, _f0JPC, _etaJPC, _JPCLS_EtaTof0Eta);
  fillJPCLS(_etaJPC, _f2JPC, _etaJPC, _JPCLS_EtaTof2Eta);
  fillJPCLS(_eta2JPC, _a0JPC, _piJPC, _JPCLS_Eta2Toa0Pi);
  fillJPCLS(_f1JPC, _a0JPC, _piJPC,_JPCLS_F1Toa0Pi);
  fillJPCLS(_etaJPC, _a2JPC, _piJPC,_JPCLS_EtaToa2Pi);
  fillJPCLS(_eta2JPC, _a2JPC, _piJPC,_JPCLS_Eta2Toa2Pi);
  fillJPCLS(_eta2JPC, _f0JPC, _etaJPC,_JPCLS_Eta2Tof0Eta);
  fillJPCLS(_eta2JPC, _f2JPC, _etaJPC,_JPCLS_Eta2Tof2Eta);
  fillJPCLS(_f1JPC, _a2JPC, _piJPC, _JPCLS_F1Toa2Pi);
  fillJPCLS(_f1JPC, _f0JPC, _etaJPC,  _JPCLS_F1Tof0Eta);
  fillJPCLS(_f1JPC, _f2JPC, _etaJPC,  _JPCLS_F1Tof2Eta);
}

JpsiGamEtaPiPiStates::~JpsiGamEtaPiPiStates()
{
}



void JpsiGamEtaPiPiStates::print(std::ostream& os) const
{
  JpsiToXGamStates::print(os);
  os << "*** eta -> a0 pi:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_EtaToa0Pi );

  os << "*** eta -> f0 eta:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_EtaTof0Eta );

  os << "*** eta -> f2 eta:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_EtaTof2Eta );

  os << "*** eta2 -> a0 pi:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_Eta2Toa0Pi );

  os << "*** eta -> a2 pi:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_EtaToa2Pi );

  os << "*** eta2 -> a2 pi:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_Eta2Toa2Pi );

  os << "*** eta2 -> f0 eta:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_Eta2Tof0Eta );

  os << "*** eta2 -> f2 eta:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_Eta2Tof2Eta );

  os << "*** f1 -> a0 pi:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_F1Toa0Pi );

  os << "*** f1 -> a2 pi:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_F1Toa2Pi );

  os << "*** f1 -> f0 eta:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_F1Tof0Eta );

  os << "*** f1 -> f2 eta:  LS combinations for the decay *** "<< std::endl;
  printDecayJPCLS(os, _JPCLS_F1Tof2Eta );  
}


