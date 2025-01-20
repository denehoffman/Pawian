//************************************************************************//
//									  //
//  Copyright 2024 Bertram Kopf (bertram@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
//									  //
//  This file is part of Pawian.					  //
//									  //
//  Pawian is free software: you can redistribute it and/or modify	  //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	  //
//									  //
//  Pawian is distributed in the hope that it will be useful,		  //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  //
//  GNU General Public License for more details.	      		  //
//									  //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	  //
//									  //
//************************************************************************//

// Fixed parametrization for the sigma pole based on
// Phys. Lett B 598 (204) 149-158

#include "PwaDynamics/SigmaParameterization.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Utils/PawianConstants.hh"
#include "qft++Extension/PawianUtils.hh"

//#include <complex>
//using std::complex;

SigmaParameterization::SigmaParameterization() :
  _mPole(0.9264)
  ,_g4pi(0.0024)
  ,_b1(0.5843)
  ,_b2(1.6663)
  ,_a(1.082)
  ,_sAdler(PawianConstants::mPiSq/2.)
  ,_mPoleSqr(0.9264*0.9264)
  ,_g2piDenom(0.9264*0.9264-PawianConstants::mPiSq/2.)
  ,_rho2pi_mPole(PawianQFT::phaseSpaceFacDefault(0.9264,PawianConstants::mPi, PawianConstants::mPi))
  ,_rho4pi_mPole(php4pi(0.9264*0.9264))
{
}
SigmaParameterization::SigmaParameterization(double mPole, double g4pi, double b1, double b2, double a) :
  _mPole(mPole)
  ,_g4pi(g4pi)
  ,_b1(b1)
  ,_b2(b2)
  ,_a(a)
  ,_sAdler(PawianConstants::mPiSq/2.)
  ,_mPoleSqr(mPole*mPole)
  ,_g2piDenom(mPole*mPole-PawianConstants::mPiSq/2.)
  ,_rho2pi_mPole(PawianQFT::phaseSpaceFacDefault(mPole,PawianConstants::mPi, PawianConstants::mPi))
  ,_rho4pi_mPole(php4pi(mPole*mPole))
{
}  
SigmaParameterization::~SigmaParameterization(){
}

complex<double> SigmaParameterization::calc(double currentMass){
  //complex<double> result(1.,0.);
  double currentM2=currentMass*currentMass;
  double fs=_b1+_b2*currentM2;
  double g2pi=fs*(currentM2-_sAdler)/_g2piDenom*exp(-(currentM2-_mPoleSqr)/_a);
  std::complex<double> rho2pi=PawianQFT::phaseSpaceFacDefault(currentMass,PawianConstants::mPi, PawianConstants::mPi);
  std::complex<double> rho4pi=php4pi(currentM2);
  std::complex<double> Gamma_tot=g2pi*rho2pi/_rho2pi_mPole+_g4pi*rho4pi/_rho4pi_mPole;
  
  complex<double> result=1./(std::complex<double>(_mPoleSqr,0.)-std::complex<double>(currentM2,0.)-PawianConstants::i*_mPole*Gamma_tot);
  return result;
}

complex<double> SigmaParameterization::calc(double currentMass, double gSigma){
  return gSigma*calc(currentMass);
}
  std::complex<double> SigmaParameterization::php4pi(double currentMass2){
  // if(currentMass2<16.*PawianConstants::mPiSq){
  // Alert << "current mass^2 must be larger than 16*m_pi^2!!!\n current mass^2: " << currentMass2 <<"\n16*m_pi^2: " << 16.*PawianConstants::mPiSq << endmsg;
  // exit(1);
  // }
  std::complex<double> result(0.,0.);
  if(currentMass2<16.*PawianConstants::mPiSq) return result;
  result=std::complex<double>(sqrt(1.-16.*PawianConstants::mPiSq/currentMass2)/(1.+exp(2.8-currentMass2)/3.5),0.);
  return result; 
}

complex<double> SigmaParameterization::calcT(double currentMass){
  double currentM2=currentMass*currentMass;
  double fs=_b1+_b2*currentM2;
  double g2pi=fs*(currentM2-_sAdler)/_g2piDenom*exp(-(currentM2-_mPoleSqr)/_a);
  std::complex<double> rho2pi=PawianQFT::phaseSpaceFacDefault(currentMass,PawianConstants::mPi, PawianConstants::mPi);
  std::complex<double> rho4pi=php4pi(currentM2);
  std::complex<double> Gamma_tot=g2pi*rho2pi/_rho2pi_mPole+_g4pi*rho4pi/_rho4pi_mPole;
  
  complex<double> result=_mPole*Gamma_tot/(std::complex<double>(_mPoleSqr,0.)-std::complex<double>(currentM2,0.)-PawianConstants::i*_mPole*Gamma_tot);
  return result;
}

