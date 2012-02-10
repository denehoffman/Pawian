#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/FitParams2K2PiGam.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamStates.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

FitParams2K2PiGam::FitParams2K2PiGam() 
{
  Psi2STo2K2PiGamStates theStates;
  theStates.print(std::cout);

  _jpclsMap[paramEnum2K2PiGam::ChiGam]=theStates.PsiToChiGamStates();
  _jpclsMap[paramEnum2K2PiGam::K892K892]=theStates.ChiTo2K892States();
  _jpclsMap[paramEnum2K2PiGam::KappaKappa]=theStates.ChiTo2K_0_States();
  _jpclsMap[paramEnum2K2PiGam::K_0_1430K_0_1430]=theStates.ChiTo2K_0_States();
  _jpclsMap[paramEnum2K2PiGam::K_2_1430K_2_1430]=theStates.ChiTo2K_2_1430States();
  _jpclsMap[paramEnum2K2PiGam::K_0_1430K_2_1430]=theStates.ChiToK0K2_States();
  _jpclsMap[paramEnum2K2PiGam::ChiToK_0_1430_K892]=theStates.ChiToKst0Kst1States();
  _jpclsMap[paramEnum2K2PiGam::K_1_1410K_1_1410]=theStates.ChiTo2K892States();
  _jpclsMap[paramEnum2K2PiGam::ChiToK_1_1410_K892]=theStates.ChiTo2K892States();
  _jpclsMap[paramEnum2K2PiGam::ChiToK_1_1400K]=theStates.ChiToK1400ToK892piStates();
  _jpclsMap[paramEnum2K2PiGam::K_1_1400ToK892Pi]=theStates.K1400ToKst1PiStates();
  _jpclsMap[paramEnum2K2PiGam::ChiToK_1_1270_K]=theStates.ChiToK1400ToK892piStates();
  _jpclsMap[paramEnum2K2PiGam::K_1_1270ToK892Pi]=theStates.K1400ToKst1PiStates();
  _jpclsMap[paramEnum2K2PiGam::K_1_1270ToK_0_1430Pi]=theStates.K1ToK0PiStates();
  _jpclsMap[paramEnum2K2PiGam::ChiToK_1_1650K]=theStates.ChiToK1400ToK892piStates(); 
  _jpclsMap[paramEnum2K2PiGam::K_1_1650ToK892Pi]=theStates.K1400ToKst1PiStates();
  _jpclsMap[paramEnum2K2PiGam::K_1_1650ToK_0_1430Pi]=theStates.K1ToK0PiStates();
  _jpclsMap[paramEnum2K2PiGam::ChiToKappaK_0_1430]=theStates.ChiTo2K_0_States();
  _jpclsMap[paramEnum2K2PiGam::K_0_1430K_0_1950]=theStates.ChiTo2K_0_States();
  _jpclsMap[paramEnum2K2PiGam::KappaK_0_1950]=theStates.ChiTo2K_0_States();
  _jpclsMap[paramEnum2K2PiGam::f980_pif1710_k]=theStates.ChiTof0f0States();
  _jpclsMap[paramEnum2K2PiGam::f980_kf1710_pi]=theStates.ChiTof0f0States();
  _jpclsMap[paramEnum2K2PiGam::f980f980]=theStates.ChiTof0f0States();
  _jpclsMap[paramEnum2K2PiGam::f980f2200]=theStates.ChiTof0f0States();
  _jpclsMap[paramEnum2K2PiGam::ChiTof980f_2_2300]=theStates.ChiTof0f2States();
  _jpclsMap[paramEnum2K2PiGam::ChiTof_2_2300sigma]=theStates.ChiTof0f2States();
  _jpclsMap[paramEnum2K2PiGam::f980_pif1370_k]=theStates.ChiTof0f0States();
  _jpclsMap[paramEnum2K2PiGam::f980_kf1370_pi]=theStates.ChiTof0f0States();
  _jpclsMap[paramEnum2K2PiGam::f980_pif1500_k]=theStates.ChiTof0f0States();
  _jpclsMap[paramEnum2K2PiGam::f980_kf1500_pi]=theStates.ChiTof0f0States();
  _jpclsMap[paramEnum2K2PiGam::f1710_pif1370_k]=theStates.ChiTof0f0States();
  _jpclsMap[paramEnum2K2PiGam::f1710_kf1370_pi]=theStates.ChiTof0f0States();
  _jpclsMap[paramEnum2K2PiGam::f980_pif_2_1270_k]=theStates.ChiTof0f2States();
  _jpclsMap[paramEnum2K2PiGam::f980_kf_2_1270_pi]=theStates.ChiTof0f2States();
  _jpclsMap[paramEnum2K2PiGam::f980_pif_2_1430_k]=theStates.ChiTof0f2States();
  _jpclsMap[paramEnum2K2PiGam::f980_kf_2_1430_pi]=theStates.ChiTof0f2States();
  _jpclsMap[paramEnum2K2PiGam::f980_pif_2_1525_k]=theStates.ChiTof0f2States();
  _jpclsMap[paramEnum2K2PiGam::f980_kf_2_1525_pi]=theStates.ChiTof0f2States();
  _jpclsMap[paramEnum2K2PiGam::f980_pif_2_1950_k]=theStates.ChiTof0f2States();
  _jpclsMap[paramEnum2K2PiGam::f980_kf_2_1950_pi]=theStates.ChiTof0f2States();
  _jpclsMap[paramEnum2K2PiGam::f1500_pif_2_1525_k]=theStates.ChiTof0f2States();
  _jpclsMap[paramEnum2K2PiGam::f1500_kf_2_1525_pi]=theStates.ChiTof0f2States();
  _jpclsMap[paramEnum2K2PiGam::f1710_pif_2_1430_k]=theStates.ChiTof0f2States();
  _jpclsMap[paramEnum2K2PiGam::f1710_kf_2_1430_pi]=theStates.ChiTof0f2States();
  _jpclsMap[paramEnum2K2PiGam::f1710_pif_2_1950_k]=theStates.ChiTof0f2States();
  _jpclsMap[paramEnum2K2PiGam::f1710_kf_2_1950_pi]=theStates.ChiTof0f2States();
  _jpclsMap[paramEnum2K2PiGam::K_0_2400KToKf980]=theStates.ChiToK0K0States();
  _jpclsMap[paramEnum2K2PiGam::ChiToK_0_2400ToKf_0_1710]=theStates.ChiToK0K0States();
  _jpclsMap[paramEnum2K2PiGam::ChiToK_1_2400K]=theStates.ChiToK1400ToK892piStates();
  _jpclsMap[paramEnum2K2PiGam::K_1_2400Tof980K]=theStates.K1pTof0KStates();
  _jpclsMap[paramEnum2K2PiGam::K892K_1_1680]=theStates.ChiTo2K892States();
  _jpclsMap[paramEnum2K2PiGam::ChiToK1680K1680]=theStates.ChiTo2K892States();
  _jpclsMap[paramEnum2K2PiGam::ChiToK1680K_0_1430]=theStates.ChiToKst0Kst1States();
  _jpclsMap[paramEnum2K2PiGam::K892K_1_2300]=theStates.ChiTo2K892States();
  _jpclsMap[paramEnum2K2PiGam::sigmaf980]=theStates.ChiTof0f0States();
  _jpclsMap[paramEnum2K2PiGam::ChiToSigmaf1370]=theStates.ChiTof0f0States();
  _jpclsMap[paramEnum2K2PiGam::sigmaf1710]=theStates.ChiTof0f0States();
  _jpclsMap[paramEnum2K2PiGam::sigmaf2200]=theStates.ChiTof0f0States();
  _jpclsMap[paramEnum2K2PiGam::K_0_1460ToK892Pi]=theStates.ChiToK0K0States();
  _jpclsMap[paramEnum2K2PiGam::K_0_1460ToK_0_1430Pi]=theStates.ChiToK0K0States();
  _jpclsMap[paramEnum2K2PiGam::K_0_1830ToK892Pi]=theStates.ChiToK0K0States();
  _jpclsMap[paramEnum2K2PiGam::K_0_1830ToK_0_1430Pi]=theStates.ChiToK0K0States();
  _jpclsMap[paramEnum2K2PiGam::ChiToPi_2_1670Pi]=theStates.ChiToPi_2PiStates();
  _jpclsMap[paramEnum2K2PiGam::Pi_2_1670Tof_2_1270Pi]=theStates.Pi_2Tof_2PiStates();
  _jpclsMap[paramEnum2K2PiGam::Pi_2_1670ToK892K]=theStates.Pi_2ToKst1KStates();
  _jpclsMap[paramEnum2K2PiGam::ChiToPi1800Pi0Tof980]=theStates.ChiToPi0Pi0States();
  _jpclsMap[paramEnum2K2PiGam::ChiToPi1800Pi0Tof1370]=theStates.ChiToPi0Pi0States();
  _jpclsMap[paramEnum2K2PiGam::ChiToPi1800Pi0ToKappa]=theStates.ChiToPi0Pi0States();
  _jpclsMap[paramEnum2K2PiGam::ChiToPi1800Pi0ToK892K]=theStates.Pi0pToKstKStates();
  _jpclsMap[paramEnum2K2PiGam::ChiToPi3000Pi0ToK892K]=theStates.Pi0pToKstKStates();
  _jpclsMap[paramEnum2K2PiGam::ChiToPi3000Pi0ToK_0_1950K]=theStates.ChiToPi0Pi0States();
  _jpclsMap[paramEnum2K2PiGam::ChiToPi_2_2285Pi]=theStates.ChiToPi_2PiStates();
  _jpclsMap[paramEnum2K2PiGam::Pi_2_2285Tof1700Pi]=theStates.Pi_2Tof0PiStates();
  _jpclsMap[paramEnum2K2PiGam::Pi_2_2285ToK892KPi]=theStates.Pi_2ToKst1KStates();
  _jpclsMap[paramEnum2K2PiGam::Pi_2_2285ToK_0_1430K]=theStates.Pi_2ToKst0KStates();
  _jpclsMap[paramEnum2K2PiGam::Pi_2_2285ToK_2_1430K]=theStates.Pi_2ToKst2KStates();
  _jpclsMap[paramEnum2K2PiGam::ChiToK_2_1770K]=theStates.ChiToK2mK0mStates();
  _jpclsMap[paramEnum2K2PiGam::K_2_1770ToK_2_1430Pi]=theStates.K2mToK2pPiStates();
  _jpclsMap[paramEnum2K2PiGam::ChiToK_0_1430KPi]=theStates.ChiTo2K_0_States();
  _jpclsMap[paramEnum2K2PiGam::ChiToK892KPi]=theStates.ChiToK1400ToK892piStates();
}

FitParams2K2PiGam::~FitParams2K2PiGam()
{
}

std::vector< boost::shared_ptr<const JPCLS> > FitParams2K2PiGam::jpclsVec(unsigned int index){
  unsigned int nOfAmps=paramEnum2K2PiGam::nAmps;
  if ( index > paramEnum2K2PiGam::nAmps ) Alert << "index " << index 
						 << "  number of amplitude parametes paramEnum2K2PiGam::numAmps = " 
 						    << nOfAmps << endmsg;

  return _jpclsMap[index]; 

}

std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& FitParams2K2PiGam::ampMap(param2K2PiGam& params, unsigned int index){
  unsigned int nOfAmps=paramEnum2K2PiGam::nAmps;
  if ( index > paramEnum2K2PiGam::nAmps ) Alert << "index " << index 
						 << "  number of amplitude parametes paramEnum2K2PiGam::numAmps = " 
 						    << nOfAmps << endmsg;

  if (index==paramEnum2K2PiGam::ChiGam) return params.PsiToChiGam;
  else if (index==paramEnum2K2PiGam::K892K892) return params.ChiTo2K892;
  else if (index==paramEnum2K2PiGam::KappaKappa) return params.ChiTo2Kappa;
  else if (index==paramEnum2K2PiGam::K_0_1430K_0_1430) return params.ChiTo2K_0_1430;
  else if (index==paramEnum2K2PiGam::K_2_1430K_2_1430) return params.ChiTo2K_2_1430;
  else if (index==paramEnum2K2PiGam::K_0_1430K_2_1430) return params.ChiToK_0_1430_K_2_1430;
  else if (index==paramEnum2K2PiGam::ChiToK_0_1430_K892) return params.ChiToK_0_1430_K892;
  else if (index==paramEnum2K2PiGam::K_1_1410K_1_1410) return params.ChiToK_1_1410_K_1_1410;
  else if (index==paramEnum2K2PiGam::ChiToK_1_1410_K892) return params.ChiToK_1_1410_K892;
  else if (index==paramEnum2K2PiGam::ChiToK_1_1400K) return params.ChiToK_1_1400K;
  else if (index==paramEnum2K2PiGam::K_1_1400ToK892Pi) return params.K1400ToK892Pi;
  else if (index==paramEnum2K2PiGam::ChiToK_1_1270_K) return params.ChiToK_1_1270_K;
  else if (index==paramEnum2K2PiGam::K_1_1270ToK892Pi) return params.K_1_1270ToK892Pi;
  else if (index==paramEnum2K2PiGam::K_1_1270ToK_0_1430Pi) return params.K_1_1270ToK_0_1430Pi;
  else if (index==paramEnum2K2PiGam::ChiToK_1_1650K) return params.ChiToK_1_1650K;
  else if (index==paramEnum2K2PiGam::K_1_1650ToK892Pi) return params.K_1_1650ToK892Pi;
  else if (index==paramEnum2K2PiGam::K_1_1650ToK_0_1430Pi) return params.K_1_1650ToK_0_1430Pi;
  else if (index==paramEnum2K2PiGam::ChiToKappaK_0_1430) return params.ChiToKappaK_0_1430;
  else if (index==paramEnum2K2PiGam::K_0_1430K_0_1950) return params.ChiToK_0_1430K_0_1950; 
  else if (index==paramEnum2K2PiGam::KappaK_0_1950) return params.ChiToKappaK_0_1950;  
  else if (index==paramEnum2K2PiGam::f980_pif1710_k) return params.f980_pif1710_k;
  else if (index==paramEnum2K2PiGam::f980_kf1710_pi) return params.f980_kf1710_pi;
  else if (index==paramEnum2K2PiGam::f980f980) return params.ChiTof980f980;
  else if (index==paramEnum2K2PiGam::f980f2200) return params.ChiTof980f2200;
  else if (index==paramEnum2K2PiGam::ChiTof980f_2_2300) return params.ChiTof980f_2_2300;
  else if (index==paramEnum2K2PiGam::ChiTof_2_2300sigma) return params.ChiTof_2_2300sigma;
  else if (index==paramEnum2K2PiGam::f980_pif1370_k) return params.f980_pif1370_k;
  else if (index==paramEnum2K2PiGam::f980_kf1370_pi) return params.f980_kf1370_pi;
  else if (index==paramEnum2K2PiGam::f980_pif1500_k) return params.f980_pif1500_k;
  else if (index==paramEnum2K2PiGam::f980_kf1500_pi) return params.f980_kf1500_pi;
  else if (index==paramEnum2K2PiGam::f1710_pif1370_k) return params.f1710_pif1370_k;
  else if (index==paramEnum2K2PiGam::f1710_kf1370_pi) return params.f1710_kf1370_pi;
  else if (index==paramEnum2K2PiGam::f980_pif_2_1270_k) return params.f980_pif_2_1270_k;
  else if (index==paramEnum2K2PiGam::f980_kf_2_1270_pi) return params.f980_kf_2_1270_pi;
  else if (index==paramEnum2K2PiGam::f980_pif_2_1430_k) return params.f980_pif_2_1430_k;
  else if (index==paramEnum2K2PiGam::f980_kf_2_1430_pi) return params.f980_kf_2_1430_pi;
  else if (index==paramEnum2K2PiGam::f980_pif_2_1525_k) return params.f980_pif_2_1525_k;
  else if (index==paramEnum2K2PiGam::f980_kf_2_1525_pi) return params.f980_kf_2_1525_pi;
  else if (index==paramEnum2K2PiGam::f980_pif_2_1950_k) return params.f980_pif_2_1950_k;
  else if (index==paramEnum2K2PiGam::f980_kf_2_1950_pi) return params.f980_kf_2_1950_pi;
  else if (index==paramEnum2K2PiGam::f1500_pif_2_1525_k) return params.f1500_pif_2_1525_k;
  else if (index==paramEnum2K2PiGam::f1500_kf_2_1525_pi) return params.f1500_kf_2_1525_pi;
  else if (index==paramEnum2K2PiGam::f1710_pif_2_1430_k) return params.f1710_pif_2_1430_k;
  else if (index==paramEnum2K2PiGam::f1710_kf_2_1430_pi) return params.f1710_kf_2_1430_pi;
  else if (index==paramEnum2K2PiGam::f1710_pif_2_1950_k) return params.f1710_pif_2_1950_k;
  else if (index==paramEnum2K2PiGam::f1710_kf_2_1950_pi) return params.f1710_kf_2_1950_pi;
  else if (index==paramEnum2K2PiGam::K_0_2400KToKf980) return params.ChiToK_0_2400ToKf980;
  else if (index==paramEnum2K2PiGam::ChiToK_0_2400ToKf_0_1710) return params.ChiToK_0_2400ToKf_0_1710;
  else if (index==paramEnum2K2PiGam::ChiToK_1_2400K) return params.ChiToK_1_2400K;
  else if (index==paramEnum2K2PiGam::K_1_2400Tof980K) return params.K_1_2400Tof980K;
  else if (index==paramEnum2K2PiGam::K892K_1_1680) return params.ChiToK892K1680;
  else if (index==paramEnum2K2PiGam::ChiToK1680K1680) return params.ChiToK1680K1680;
  else if (index==paramEnum2K2PiGam::ChiToK1680K_0_1430) return params.ChiToK1680K_0_1430;
  else if (index==paramEnum2K2PiGam::K892K_1_2300) return params.ChiToK892K2300;
  else if (index==paramEnum2K2PiGam::sigmaf980) return params.ChiToSigmaf980;
  else if (index==paramEnum2K2PiGam::ChiToSigmaf1370) return params.ChiToSigmaf1370;
  else if (index==paramEnum2K2PiGam::sigmaf1710) return params.ChiToSigmaf1710;
  else if (index==paramEnum2K2PiGam::sigmaf2200) return params.ChiToSigmaf2200;
  else if (index==paramEnum2K2PiGam::K_0_1460ToK892Pi) return params.K_0_1460ToK892Pi;
  else if (index==paramEnum2K2PiGam::K_0_1460ToK_0_1430Pi) return params.K_0_1460ToK_0_1430Pi;
  else if (index==paramEnum2K2PiGam::K_0_1830ToK892Pi) return params.K_0_1830ToK892Pi;
  else if (index==paramEnum2K2PiGam::K_0_1830ToK_0_1430Pi) return params.K_0_1830ToK_0_1430Pi;
  else if (index==paramEnum2K2PiGam::ChiToPi_2_1670Pi) return params.ChiToPi_2_1670Pi;
  else if (index==paramEnum2K2PiGam::Pi_2_1670Tof_2_1270Pi) return params.Pi_2_1670Tof_2_1270Pi;
  else if (index==paramEnum2K2PiGam::Pi_2_1670ToK892K) return params.Pi_2_1670ToK892K;
  else if (index==paramEnum2K2PiGam::ChiToPi1800Pi0Tof980) return params.ChiToPi1800Pi0Tof980;
  else if (index==paramEnum2K2PiGam::ChiToPi1800Pi0Tof1370) return params.ChiToPi1800Pi0Tof1370;
  else if (index==paramEnum2K2PiGam::ChiToPi1800Pi0ToKappa) return params.ChiToPi1800Pi0ToKappa;
  else if (index==paramEnum2K2PiGam::ChiToPi1800Pi0ToK892K) return params.ChiToPi1800Pi0ToK892K;
  else if (index==paramEnum2K2PiGam::ChiToPi3000Pi0ToK892K) return params.ChiToPi3000Pi0ToK892K;
  else if (index==paramEnum2K2PiGam::ChiToPi3000Pi0ToK_0_1950K) return params.ChiToPi3000Pi0ToK_0_1950K;
  else if (index==paramEnum2K2PiGam::ChiToPi_2_2285Pi) return params.ChiToPi_2_2285Pi; 
  else if (index==paramEnum2K2PiGam::Pi_2_2285Tof1700Pi) return params.Pi_2_2285Tof1700Pi;
  else if (index==paramEnum2K2PiGam::Pi_2_2285ToK892KPi) return params.Pi_2_2285ToK892KPi;
  else if (index==paramEnum2K2PiGam::Pi_2_2285ToK_0_1430K) return params.Pi_2_2285ToK_0_1430K;
  else if (index==paramEnum2K2PiGam::Pi_2_2285ToK_2_1430K) return params.Pi_2_2285ToK_2_1430K;
  else if (index==paramEnum2K2PiGam::ChiToK_2_1770K) return params.ChiToK_2_1770K;
  else if (index==paramEnum2K2PiGam::K_2_1770ToK_2_1430Pi) return params.K_2_1770ToK_2_1430Pi;
  else if (index==paramEnum2K2PiGam::ChiToK_0_1430KPi) return params.ChiToK_0_1430KPi;
  else if (index==paramEnum2K2PiGam::ChiToK892KPi) return params.ChiToK892KPi; 
  else{
    Alert << "index " << index << " not supported !!!" << endmsg;
    exit(0);
  }

  return params.PsiToChiGam; //makes the compiler happy
}

pair<double, double>& FitParams2K2PiGam::massPair(param2K2PiGam& params, unsigned int index){
  unsigned int nOfAmps=paramEnum2K2PiGam::nAmps;
  unsigned int nOfMasses=paramEnum2K2PiGam::nMasses;
  if ( index<paramEnum2K2PiGam::nAmps || index > paramEnum2K2PiGam::nMasses ) Alert << "index " << index 
										      << "  must be between  " << nOfAmps
										      << " and " << nOfMasses << endmsg;
  if (index==paramEnum2K2PiGam::K892) return params.BwK892;
  else if (index==paramEnum2K2PiGam::Kappa) return params.BwKappa;
  else if (index==paramEnum2K2PiGam::K_0_1430) return params.BwK_0_1430;
  else if (index==paramEnum2K2PiGam::K_1_1400) return params.BwK_1_1400;
  else if (index==paramEnum2K2PiGam::K_1_1410) return params.BwK_1_1410;
  else if (index==paramEnum2K2PiGam::K_2_1430) return params.BwK_2_1430;
  else if (index==paramEnum2K2PiGam::K_1_1270) return params.BwK_1_1270;
  else if (index==paramEnum2K2PiGam::K_1_1650) return params.BwK_1_1650;
  else if (index==paramEnum2K2PiGam::f1500) return params.Bwf1500;
  else if (index==paramEnum2K2PiGam::f1710) return params.Bwf1710;
  else if (index==paramEnum2K2PiGam::f2200) return params.Bwf2200;
  else if (index==paramEnum2K2PiGam::sigma) return params.BwSigma;
  else if (index==paramEnum2K2PiGam::f1360) return params.Bwf1360;
  else if (index==paramEnum2K2PiGam::f1370) return params.Bwf1370;
  else if (index==paramEnum2K2PiGam::f_2_1270) return params.Bwf_2_1270;
  else if (index==paramEnum2K2PiGam::f_2_1430) return params.Bwf_2_1430;
  else if (index==paramEnum2K2PiGam::f_2_1525) return params.Bwf_2_1525;
  else if (index==paramEnum2K2PiGam::f_2_1950) return params.Bwf_2_1950;
  else if (index==paramEnum2K2PiGam::f_2_2300) return params.Bwf_2_2300;
  else if (index==paramEnum2K2PiGam::K_0_2400) return params.BwK_0_2400;
  else if (index==paramEnum2K2PiGam::K_1_2400) return params.BwK_1_2400;
  else if (index==paramEnum2K2PiGam::K_0_1950) return params.BwK_0_1950;
  else if (index==paramEnum2K2PiGam::K_1_1680) return params.BwK_1_1680;
  else if (index==paramEnum2K2PiGam::K_1_2300) return params.BwK_1_2300;
  else if (index==paramEnum2K2PiGam::K_0_1460) return params.BwK_0_1460;
  else if (index==paramEnum2K2PiGam::K_0_1830) return params.BwK_0_1830;
  else if (index==paramEnum2K2PiGam::Pi_2_1670) return params.BwPi_2_1670;
  else if (index==paramEnum2K2PiGam::Pi1800) return params.BwPi1800;
  else if (index==paramEnum2K2PiGam::Pi3000) return params.BwPi3000;
  else if (index==paramEnum2K2PiGam::Pi_2_2285) return params.BwPi_2_2285;
  else if (index==paramEnum2K2PiGam::K_2_1770) return params.BwK_2_1770;
  else{ Alert << "index " << index << " not supported !!!" << endmsg;
    exit(0);
  }

  return params.BwK_0_1430; //makes the compiler happy
}

int FitParams2K2PiGam::setFitParamValDec(param2K2PiGam& theParamVal, const std::vector<double>& par, int counter, unsigned int index){

  int resultCount=counter;
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > currentStates;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >* currentMap=0;

  currentStates = _jpclsMap[index];
  currentMap = &ampMap(theParamVal, index);

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;

  for ( itJPCLS=currentStates.begin(); itJPCLS!=currentStates.end(); ++itJPCLS){
    double mag=par[resultCount];
    resultCount++;
    
    double phi=par[resultCount];
    resultCount++;
    

    std::pair <double,double> tmpParameter=make_pair(mag,phi);
    (*currentMap)[(*itJPCLS)]=tmpParameter;
  }

  return resultCount;
}


int FitParams2K2PiGam::setFitParamValMass(param2K2PiGam& theParamVal, const std::vector<double>& par, int counter, unsigned int index){

  int resultCount=counter;

  pair<double, double>* currentPair= & massPair(theParamVal, index);


  double currentMass=par[resultCount];
  resultCount++;
  double currentWidth=par[resultCount];
  resultCount++;

  std::pair <double,double> thePair=make_pair(currentMass, currentWidth);

  (*currentPair)=thePair;

  return resultCount;
}


int FitParams2K2PiGam::setFitParamFlatteMass(param2K2PiGam& theParamVal, const std::vector<double>& par, int counter, std::string key){

  int resultCount=counter;

  if (key==paramEnum2K2PiGam::name(paramEnum2K2PiGam::f980)){ 
    double* resultFlatMf980= &theParamVal.Flatf980;
    double* resultFlatgKK= &theParamVal.Flatf980gKK;
    double* resultFlatgPiPi= &theParamVal.Flatf980gPiPi;
    
    double currentMass=par[resultCount];
    resultCount++;
    (*resultFlatMf980)=currentMass;
    
    double currentFlatgKK=par[resultCount];
    resultCount++;
    (*resultFlatgKK)=currentFlatgKK;
    
    double currentFlatgPiPi=par[resultCount];
    resultCount++;
    (*resultFlatgPiPi)=currentFlatgPiPi;
  }
  else if (key==paramEnum2K2PiGam::name(paramEnum2K2PiGam::phaseSpace)){
    double* resultPhaseSpace= &theParamVal.phaseSpace;
    double currentPhaseSpace=par[resultCount];
    resultCount++;
    (*resultPhaseSpace)=currentPhaseSpace;
  }
  else { Alert << "Key: " << key << " not supported for setting up the Flatte mass parameters!!!" << endmsg;
    exit(1);
  }


  return resultCount;
}

void FitParams2K2PiGam::setMnUsrParamsDec(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal, unsigned int index){
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > startParams;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > errParams;

  startParams= ampMap(startVal, index);
  errParams= ampMap(errVal, index);


  int counter=0;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for (it=startParams.begin(); it!=startParams.end(); ++it){
    boost::shared_ptr<const JPCLS> theJPCLS=it->first;   
    std::pair <double,double> startPair=it->second;
    std::pair <double,double> errPair=errParams.find(theJPCLS)->second;

    //now fill the fitParameterMap
    std::string magStr=theJPCLS->name()+paramEnum2K2PiGam::name(index)+"Mag";
    std::string phiStr=theJPCLS->name()+paramEnum2K2PiGam::name(index)+"phi";

    double magVal=startPair.first;
    double phiVal=startPair.second;
    
    double magErr=errPair.first;
    double phiErr=errPair.second;
    
    double magMin=magVal-magErr;
    if (magMin<0.) magMin=0.;
    
    upar.Add(magStr, magVal, magErr, magMin, magVal+magErr);
    upar.Add(phiStr, phiVal, phiErr, -3*M_PI, 3*M_PI);

    counter++;
  }

}


void FitParams2K2PiGam::setMnUsrParamsMass(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal, unsigned int index){
  pair<double, double> startParams;
  pair<double, double> errParams;
  startParams=massPair(startVal, index);
  errParams=massPair(errVal, index);

  //now fill the fitParameterMap
  std::string massStr=paramEnum2K2PiGam::name(index)+"Mag";
  std::string widthStr=paramEnum2K2PiGam::name(index)+"Width";
  
  double massVal=startParams.first;
  double widthVal=startParams.second;
  
  double massErr=errParams.first;
  double widthErr=errParams.second;

  double massMin=massVal-3.*massErr;
  if (massMin<0.) massMin=0.;

  double massMax=massVal+3.*massErr;

  double widthMin=0.; 
  double widthMax=2*widthVal;
 
  upar.Add(massStr, massVal, massErr, massMin, massMax);
  upar.Add(widthStr, widthVal, widthErr, widthMin, widthMax);
}


void FitParams2K2PiGam::setMnUsrParamsFlatteMass(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal, std::string key){
 
  if (key==paramEnum2K2PiGam::name(paramEnum2K2PiGam::f980))
    {
      double start_f980M=startVal.Flatf980;
      double start_FlatgKK=startVal.Flatf980gKK;
      double start_FlatgPiPi=startVal.Flatf980gPiPi;
      double err_f980M=errVal.Flatf980;
      double err_FlatgKK=errVal.Flatf980gKK;
      double err_FlatgPiPi=errVal.Flatf980gPiPi;
      
      //now fill the fitParameterMap
      std::string massStr=paramEnum2K2PiGam::name(paramEnum2K2PiGam::f980);
      std::string gKKStr=paramEnum2K2PiGam::name(paramEnum2K2PiGam::f980gKK);
      std::string gPiPiStr=paramEnum2K2PiGam::name(paramEnum2K2PiGam::f980gpipi);
      
      double massMin=0.97;
      double massMax=1.2;
      upar.Add(massStr, start_f980M, err_f980M, massMin, massMax);
      upar.Add(gKKStr, start_FlatgKK, err_FlatgKK, 0., start_FlatgKK+3*err_FlatgKK);
      upar.Add(gPiPiStr, start_FlatgPiPi, err_FlatgPiPi, 0., start_FlatgPiPi+3*err_FlatgPiPi);
    }
  
  else if ( key == paramEnum2K2PiGam::name(paramEnum2K2PiGam::phaseSpace)){ 
    double start_phaseSpace=startVal.phaseSpace;
    double err_phaseSpace=errVal.phaseSpace;
    std::string theStr=paramEnum2K2PiGam::name(paramEnum2K2PiGam::phaseSpace);
    
    double phaseSpaceMax=start_phaseSpace+3.*err_phaseSpace;
    if ( phaseSpaceMax > 1. ) phaseSpaceMax=1.;
    upar.Add(theStr, start_phaseSpace, err_phaseSpace, 0., phaseSpaceMax);
  }
  
  else { Alert << "Key: " << key << " not supported for setting up the Flatte mass parameters!!!" << endmsg;
    exit(1);
  } 
}
