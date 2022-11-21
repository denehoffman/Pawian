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

#include "gslUtils/CMIntegration.hh"
#include "Utils/Faddeeva.hh"
#include "qft++Extension/PawianUtils.hh"
#include "Utils/PawianConstants.hh"

#include <complex>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

CMunstable_params CMIntegration::_CMunstable_params(0.90338253, 0.59405387, 0.1349768, 0.13957, 0.13957);
bool CMIntegration::_calcRealPart=true;
std::complex<double> CMIntegration::_currentS(1.,0.);

CMIntegration::CMIntegration(double mpole, double fpole, double mu, double m1, double m2){
  //_CMunstable_params._mPole=mpole;
  gsl_set_error_handler_off();
}

CMIntegration::~CMIntegration(){
}

std::complex<double> CMIntegration::integrate(std::complex<double> s){
  _currentS=s;
  double epsabs = 1e-8;
  double epsrel = 1e-8;
  size_t limit = 1000;

  double intlowerborder=(_CMunstable_params._m1+_CMunstable_params._m2)*(_CMunstable_params._m1+_CMunstable_params._m2);
  // set gsl function
  gsl_function F;
  F.function = &FIntWrapper;

  const size_t n=1000;
  gsl_integration_workspace* wsp1=gsl_integration_workspace_alloc(n);

  _calcRealPart=true;
  double resultReal, abserrReal;
  gsl_integration_qagiu(&F, intlowerborder, epsabs, epsrel, limit, wsp1, &resultReal, &abserrReal); 

  //std::cout << "\ns: "<< s << std::endl;
  //std::cout << "result (real): " << resultReal << std::endl;
  //std::cout << "abserr (real): " << abserrReal << std::endl;

  _calcRealPart=false;
  double resultImag, abserrImag;
  gsl_integration_qagiu(&F, intlowerborder, epsabs, epsrel, limit, wsp1, &resultImag, &abserrImag);
  //std::cout << "result (imag): " << resultImag << std::endl;
  //std::cout << "abserr (imag): " << abserrImag << std::endl;

  complex<double> result(resultReal, resultImag);
  return result;
}

double CMIntegration::FIntWrapper(double x, void * params){
  double result=0.;
  if (_calcRealPart) result=(Ctilde(_currentS,x)).real();
  else result=(Ctilde(_currentS,x)).imag();
  return result; 
}

std::complex<double> CMIntegration::Ctilde(std::complex<double> s, double sprime){
     std::complex<double> result =
      -(1./PawianConstants::pi)*pow(_CMunstable_params._fPole,2.)
	 *(Sigma(sprime, _CMunstable_params._m1,_CMunstable_params._m2)).imag()
      /dsNorm(sprime, _CMunstable_params._m1, _CMunstable_params._m2, _CMunstable_params._mPole, _CMunstable_params._fPole)
       *PawianQFT::ChewMandelstamDudek(s, sqrt(sprime), _CMunstable_params._mu);

   return result;
  
}

std::complex<double> CMIntegration::Sigma(double sprime, double m1, double m2){
  complex<double> sprimecompl(sprime,0.);  
  std::complex<double> result = (sprime-(m1+m2)*(m1+m2))*PawianQFT::ChewMandelstamDudek(sprimecompl, m1, m2);
    return result;
}

double CMIntegration::dsNorm(std::complex<double> s, double m1, double m2, double mpole, double fpole){
  complex<double> mpolecomplex(mpole,0.);
  complex<double> fpolecomplex(fpole,0.);
  complex<double> m1complex(m1,0.);
  complex<double> m2complex(m2,0.);
  complex<double> ds=s-mpolecomplex*mpolecomplex
    +(fpolecomplex*fpolecomplex)*(s-(m1complex+m2complex)*(m1complex+m2complex))
    *PawianQFT::ChewMandelstamDudek(s, m1, m2);
  return norm(ds);
}

