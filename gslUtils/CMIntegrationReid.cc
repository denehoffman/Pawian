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

#include "gslUtils/CMIntegrationReid.hh"
#include "Utils/Faddeeva.hh"
#include "qft++Extension/PawianUtils.hh"
#include "Utils/PawianConstants.hh"
#include "ErrLogger/ErrLogger.hh"

#include <complex>
#include <gsl/gsl_errno.h>


CMIntegrationReid::CMIntegrationReid(double mpole, double fpole, double mu, double m1, double m2) :
  CMIntegration(mpole, fpole, mu, m1, m2){

  gsl_function Fnorm;
  Fnorm.function=FIntWrapperNormCondition;

  double result, abserr;  
  std::string fitName="normalizationConstant";
  doFit(Fnorm, result, abserr, fitName);
  InfoMsg << "normalization constant: " << result << " +- " << abserr << endmsg;
}



CMIntegrationReid::~CMIntegrationReid(){
}


void CMIntegrationReid::integrate(std::complex<double> s, std::complex<double>& result, std::complex<double>& resulterr){
  _currentS=s;
  if( _currentS.imag() > -1.e-12) _currentS=complex<double>(_currentS.real(), -1.e-12); //be shure to be in the 4th quadrant  

  // set gsl function
  gsl_function F;
  F.function = &FIntWrapper;
  
  _calcRealPart=true;
  double resultReal, abserrReal;
  std::string fitName="Reid real part";
  doFit(F, resultReal, abserrReal, fitName);

  
  _calcRealPart=false;
  double resultImag, abserrImag;
  fitName="Reid imaginary part";
  doFit(F, resultImag, abserrImag, fitName);
  
  result=complex<double>(resultReal, resultImag);
  resulterr=complex<double>(abserrReal, abserrImag);
  return;
}

double CMIntegrationReid::FIntWrapper(double x, void * params){
  double result=0.;
  if (_calcRealPart) result=(Ctilde(_currentS,x)).real();
  else result=(Ctilde(_currentS,x)).imag();
  return result; 
}

double CMIntegrationReid::FIntWrapperNormCondition(double x, void * params){
    double result=-(1./PawianConstants::pi)*dsImag(x, _CMunstable_params._m1, _CMunstable_params._m2, _CMunstable_params._mPole, _CMunstable_params._fPole)/dsNorm(x, _CMunstable_params._m1, _CMunstable_params._m2, _CMunstable_params._mPole, _CMunstable_params._fPole);
  // double result=-(1./PawianConstants::pi)
  //   *_CMunstable_params._fPole*_CMunstable_params._fPole
  //   *Sigma(x, _CMunstable_params._m1, _CMunstable_params._m2).imag()/dsNorm(x, _CMunstable_params._m1, _CMunstable_params._m2, _CMunstable_params._mPole, _CMunstable_params._fPole);
  return result;
}

std::complex<double> CMIntegrationReid::Ctilde(std::complex<double> s, double sprime){
  std::complex<double> result =
      -(1./PawianConstants::pi)*pow(_CMunstable_params._fPole,2.)
	 *(Sigma(sprime, _CMunstable_params._m1,_CMunstable_params._m2)).imag()
       /dsNorm(sprime, _CMunstable_params._m1, _CMunstable_params._m2, _CMunstable_params._mPole, _CMunstable_params._fPole)
       *PawianQFT::ChewMandelstamReid(s, sqrt(sprime), _CMunstable_params._mu);
   return result;  
}

std::complex<double> CMIntegrationReid::Sigma(double sprime, double m1, double m2){
  std::complex<double> sprimecompl(sprime,0.);
  std::complex<double> result = (sprime-(m1+m2)*(m1+m2))*PawianQFT::ChewMandelstamReid(sprimecompl, m1, m2);
    return result;
}

double CMIntegrationReid::dsNorm(double sprime, double m1, double m2, double mpole, double fpole){
  complex<double> sprimecomplex(sprime,0.);
  complex<double> mpolecomplex(mpole,0.);
  complex<double> fpolecomplex(fpole,0.);
  complex<double> m1complex(m1,0.);
  complex<double> m2complex(m2,0.);
  complex<double> ds=sprimecomplex-mpolecomplex*mpolecomplex
    +(fpolecomplex*fpolecomplex)*(sprime-(m1complex+m2complex)*(m1complex+m2complex))
     *PawianQFT::ChewMandelstamReid(sprimecomplex, m1, m2);
  return norm(ds);
}


double CMIntegrationReid::dsImag(double sprime, double m1, double m2, double mpole, double fpole){
  complex<double> sprimecomplex(sprime,0.);
  complex<double> mpolecomplex(mpole,0.);
  complex<double> fpolecomplex(fpole,0.);
  complex<double> m1complex(m1,0.);
  complex<double> m2complex(m2,0.);
  complex<double> ds=sprimecomplex-mpolecomplex*mpolecomplex
     +(fpolecomplex*fpolecomplex)*(sprime-(m1complex+m2complex)*(m1complex+m2complex))
     *PawianQFT::ChewMandelstamReid(sprimecomplex, m1, m2);
  return ds.imag();
}

