//************************************************************************//
//									  //
//  Copyright 2023 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Meike Kuessner (mkuessner@ep1.rub.de)		  //
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

#include "gslUtils/CMIntegrationDudek.hh"
#include "Utils/Faddeeva.hh"
#include "qft++Extension/PawianUtils.hh"
#include "Utils/PawianConstants.hh"
#include "ErrLogger/ErrLogger.hh"

#include <complex>
#include <gsl/gsl_errno.h>

std::complex<double> CMIntegrationDudek::_offset(0.,0.);

CMIntegrationDudek::CMIntegrationDudek(double mpole, double fpole, double mu, double m1, double m2, std::string offsetMethod) :
  CMIntegration(mpole, fpole, mu, m1, m2)
  ,_offsetMethod(offsetMethod){
  if(_offsetMethod != "Reid" && _offsetMethod!="Pi"){
    Alert << "offset method: " << _offsetMethod << " no supported!!!! Use Reid or Pi instead" << endmsg;
    exit(1);
  }
 
  gsl_function F;
  std::string fitNameOffset;
  if(_offsetMethod == "Reid"){
    F.function = &FIntWrapperOffsetDSigDdPi;
    InfoMsg << "use integration method: Reid offset" << endmsg;
    fitNameOffset="Dudek_ReidOffset";
  }
  else{ //_offsetMethod == "Pi"
    F.function = &FIntWrapperOffsetDSigDdReid;
    InfoMsg << "use integration method: pi offset" << endmsg;
    fitNameOffset="Dudek_piOffset";
  }

  double result, abserr;
  doFit(F, result, abserr, fitNameOffset);

  _offset=complex<double>(result, 0.);
  InfoMsg << "integration for pi: obtained value for offset: " << _offset << endmsg;


  gsl_function Fnorm;
  Fnorm.function=FIntWrapperNormCondition;
  
  std::string fitName="normalizationConstant";
  doFit(Fnorm, result, abserr, fitName);
  InfoMsg << "normalization constant: " << result << " +- " << abserr << endmsg;

  double superposMassSqr=2.*2.;
  double superpsoWeight=(1./PawianConstants::pi)*_CMunstable_params._fPole*_CMunstable_params._fPole*(Sigma(superposMassSqr, _CMunstable_params._m1, _CMunstable_params._m2)).imag()/dsNorm(superposMassSqr, _CMunstable_params._m1, _CMunstable_params._m2, _CMunstable_params._mPole, _CMunstable_params._fPole);
  InfoMsg << "superposition weight at 2 MeV: " << superpsoWeight  << endmsg;
}



CMIntegrationDudek::~CMIntegrationDudek(){
}


void CMIntegrationDudek::integrate(std::complex<double> s, std::complex<double>& result, std::complex<double>& resulterr){
  _currentS=s;
  if( _currentS.imag() > -1.e-12) _currentS=complex<double>(_currentS.real(), -1.e-12); //be shure to be in the 4th quadrant  

  // set gsl function
  gsl_function F;
  if(_offsetMethod == "Reid") F.function = &FIntWrapperPi;
  else F.function = &FIntWrapperReid;
  
  _calcRealPart=true;
  double resultReal, abserrReal;
  std::string fitName="Dudek real part";
  doFit(F, resultReal, abserrReal, fitName);

  
  _calcRealPart=false;
  double resultImag, abserrImag;
  fitName="Dudek imaginary part";
  doFit(F, resultImag, abserrImag, fitName);
  
  result=complex<double>(resultReal, resultImag)+_offset;
  resulterr=complex<double>(abserrReal, abserrImag);
  return;
}

double CMIntegrationDudek::FIntWrapperReid(double x, void * params){
  double result=0.;
  if (_calcRealPart) result=(CtildeReid(_currentS,x)).real();
  else result=(CtildeReid(_currentS,x)).imag();
  return result; 
}

double CMIntegrationDudek::FIntWrapperPi(double x, void * params){
  double result=0.;
  if (_calcRealPart) result=(CtildePi(_currentS,x)).real();
  else result=(CtildePi(_currentS,x)).imag();
  return result; 
}

double CMIntegrationDudek::FIntWrapperOffsetDSigDdPi(double x, void * params){
  double result= OffsetDSigDdPi(x);
  return result; 
}

double CMIntegrationDudek::FIntWrapperOffsetDSigDdReid(double x, void * params){
  double result = OffsetDSigDdReid(x);
  return result; 
}

double CMIntegrationDudek::FIntWrapperNormCondition(double x, void * params){
   double result=-(1./PawianConstants::pi)*dsImag(x, _CMunstable_params._m1, _CMunstable_params._m2, _CMunstable_params._mPole, _CMunstable_params._fPole)/dsNorm(x, _CMunstable_params._m1, _CMunstable_params._m2, _CMunstable_params._mPole, _CMunstable_params._fPole);
  // double result=-(1./PawianConstants::pi)
  //   *_CMunstable_params._fPole*_CMunstable_params._fPole
  //   *Sigma(x, _CMunstable_params._m1, _CMunstable_params._m2).imag()/dsNorm(x, _CMunstable_params._m1, _CMunstable_params._m2, _CMunstable_params._mPole, _CMunstable_params._fPole);
  return result;
}

std::complex<double> CMIntegrationDudek::CtildeReid(std::complex<double> s, double sprime){
  std::complex<double> result =
      (1./PawianConstants::pi)*pow(_CMunstable_params._fPole,2.)
	 *(Sigma(sprime, _CMunstable_params._m1,_CMunstable_params._m2)).imag()
       /dsNorm(sprime, _CMunstable_params._m1, _CMunstable_params._m2, _CMunstable_params._mPole, _CMunstable_params._fPole)
       *PawianQFT::ChewMandelstamReid(s, sqrt(sprime), _CMunstable_params._mu);
   return result;
  
}

std::complex<double> CMIntegrationDudek::CtildePi(std::complex<double> s, double sprime){
  complex<double> cm_s_sprime_mu= complex<double>(-1./PawianConstants::pi,0.)+PawianQFT::ChewMandelstamDudek(s, sqrt(sprime), _CMunstable_params._mu);
     std::complex<double> result =
      -(1./PawianConstants::pi)*pow(_CMunstable_params._fPole,2.)
	 *(Sigma(sprime, _CMunstable_params._m1,_CMunstable_params._m2)).imag()
       /dsNorm(sprime, _CMunstable_params._m1, _CMunstable_params._m2, _CMunstable_params._mPole, _CMunstable_params._fPole)
       *cm_s_sprime_mu;
   return result;
}


double CMIntegrationDudek::OffsetDSigDdPi(double sprime){
double result = -1./PawianConstants::pi * 1./PawianConstants::pi
  * pow(_CMunstable_params._fPole,2.)
  *(Sigma(sprime, _CMunstable_params._m1,_CMunstable_params._m2)).imag()
  /dsNorm(sprime, _CMunstable_params._m1, _CMunstable_params._m2, _CMunstable_params._mPole, _CMunstable_params._fPole);
  return result;
}

double CMIntegrationDudek::OffsetDSigDdReid(double sprime){
double result = -1./PawianConstants::pi * 1./PawianConstants::pi
  * pow(_CMunstable_params._fPole,2.)
  *(Sigma(sprime, _CMunstable_params._m1,_CMunstable_params._m2)).imag()
  /dsNorm(sprime, _CMunstable_params._m1, _CMunstable_params._m2, _CMunstable_params._mPole, _CMunstable_params._fPole)
 *(1.+(2.*sqrt(sprime)*_CMunstable_params._mu*log(sqrt(sprime)/_CMunstable_params._mu))/(sprime-_CMunstable_params._mu*_CMunstable_params._mu));
  return result;
}


std::complex<double> CMIntegrationDudek::Sigma(double sprime, double m1, double m2){
  std::complex<double> sprimecompl(sprime,0.);
  std::complex<double> result = (sprime-(m1+m2)*(m1+m2))
    *PawianQFT::ChewMandelstamDudek(sprimecompl, m1, m2);
    return result;
}

double CMIntegrationDudek::dsNorm(double sprime, double m1, double m2, double mpole, double fpole){
  complex<double> sprimecomplex(sprime,0.);
  complex<double> mpolecomplex(mpole,0.);
  complex<double> fpolecomplex(fpole,0.);
  complex<double> m1complex(m1,0.);
  complex<double> m2complex(m2,0.);
  complex<double> ds=sprimecomplex-mpolecomplex*mpolecomplex
     +(fpolecomplex*fpolecomplex)*(sprime-(m1complex+m2complex)*(m1complex+m2complex))
     *PawianQFT::ChewMandelstamDudek(sprimecomplex, m1, m2);
  return norm(ds);
}

double CMIntegrationDudek::dsImag(double sprime, double m1, double m2, double mpole, double fpole){
  complex<double> sprimecomplex(sprime,0.);
  complex<double> mpolecomplex(mpole,0.);
  complex<double> fpolecomplex(fpole,0.);
  complex<double> m1complex(m1,0.);
  complex<double> m2complex(m2,0.);
  complex<double> ds=sprimecomplex-mpolecomplex*mpolecomplex
     +(fpolecomplex*fpolecomplex)*(sprime-(m1complex+m2complex)*(m1complex+m2complex))
     *PawianQFT::ChewMandelstamDudek(sprimecomplex, m1, m2);
  return ds.imag();
}

