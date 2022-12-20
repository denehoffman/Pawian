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

#include "gslUtils/CMIntegration2DReid.hh"
#include "Utils/Faddeeva.hh"
#include "qft++Extension/PawianUtils.hh"
#include "Utils/PawianConstants.hh"
#include "ErrLogger/ErrLogger.hh"

#include <complex>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>


CMIntegration2DReid::CMIntegration2DReid(double mpole1, double fpole1, double m11, double m21,double mpole2, double fpole2, double m12, double m22) :
  CMIntegration2D(mpole1, fpole1, m11, m21, mpole2, fpole2, m12, m22)
  {
}

CMIntegration2DReid::~CMIntegration2DReid(){
}

void CMIntegration2DReid::setup(){
}

  void CMIntegration2DReid::integrate(std::complex<double> s, std::complex<double>& result, std::complex<double>& resulterr){
    _currentS=s;
    if( _currentS.imag() > -1.e-12) _currentS=complex<double>(_currentS.real(), -1.e-12); //be shure to be in the 4th quadrant

    // set gsl function
    gsl_function F;
    F.function = &FouterIntWrapper;

    _calcRealPart=true;
    double resultReal, abserrReal;
    std::string fitName="Reid real part";
    doFit2(F, resultReal, abserrReal, fitName);
    result=resultReal;
    //result=_currentOuterResult;
    
    _calcRealPart=false;
    double resultImag, abserrImag;
    fitName="Reid imaginary part";
    doFit2(F, resultImag, abserrImag, fitName);
    
    result=complex<double>(resultReal, resultImag);
    resulterr=complex<double>(abserrReal, abserrImag);
    return;
}

double CMIntegration2DReid::FouterIntWrapper(double x, void * params){
  double result=0.;
  _currentsprimeprime=x;
  //std::cout << "sprimeprime outer: " << x << std::endl;
  std::complex<double> outerResult =(1./PawianConstants::pi)*(1./PawianConstants::pi)*pow(_CMunstable_params2._fPole,2.)
    *(Sigma(_currentsprimeprime, _CMunstable_params2._m1, _CMunstable_params2._m2)).imag()
    / dsNorm(_currentsprimeprime, _CMunstable_params2._m1, _CMunstable_params2._m2, _CMunstable_params2._mPole, _CMunstable_params2._fPole);

  // set gsl function
  gsl_function F;
  F.function = &FinnerIntWrapper;  
  const size_t n=1000;
  gsl_integration_workspace* wsp2=gsl_integration_workspace_alloc(n);
  double epsabs=1.e-7;
  double epsrel=1.e-7;
  size_t limit=1000;

  double resultReal1;
  double resultReal1Err;
  _calcRealPartInner=true;
  int statusInner = gsl_integration_qagiu(&F, _integLowerBorder, epsabs, epsrel, limit, wsp2, &resultReal1, &resultReal1Err);
  gsl_integration_workspace_free (wsp2);

  double resultImag1;
  double resultImag1Err;
  _calcRealPartInner=false;
  gsl_integration_workspace* wsp2Imag=gsl_integration_workspace_alloc(n);
  statusInner = gsl_integration_qagiu(&F, _integLowerBorder, epsabs, epsrel, limit, wsp2Imag, &resultImag1, &resultImag1Err);
  gsl_integration_workspace_free (wsp2Imag);
  
  //std::cout << "outerResult: " << outerResult << std::endl;
  //std::cout << "result1: " << result1 << std::endl;
  //std::cout << "_currentInnerResult: " << _currentInnerResult << std::endl;
  //std::complex<double> resCompl=outerResult*_currentInnerResult;
  std::complex<double> result1Compl(resultReal1, resultImag1);
  std::complex<double> resCompl=outerResult*result1Compl;
  //_currentOuterResult=resCompl;

  if(_calcRealPart) result = resCompl.real();
  else result = resCompl.imag(); 
  return result; 
}

double CMIntegration2DReid::FinnerIntWrapper(double x, void * params){
  //  double *theParams = (double *) params;
  // std::cout << "_currentsprimeprime: " << _currentsprimeprime << std::endl;
  // std::cout << "sprime: " << x << std::endl;
  double result=0.;
  double sprime=x;
  std::complex<double> innerResult = pow(_CMunstable_params._fPole,2.)
    *(Sigma(sprime, _CMunstable_params._m1, _CMunstable_params._m2)).imag()
    / dsNorm(sprime, _CMunstable_params._m1, _CMunstable_params._m2, _CMunstable_params._mPole, _CMunstable_params._fPole)
    *PawianQFT::ChewMandelstamReid(_currentS, sqrt(sprime), sqrt(_currentsprimeprime));
  //std::cout << "innerResult: " << innerResult << std::endl;
  //_currentInnerResult=innerResult;

  if(_calcRealPartInner) result=innerResult.real();
  else result=innerResult.imag();
  return result;
}

std::complex<double> CMIntegration2DReid::Sigma(double sprime, double m1, double m2){
  std::complex<double> sprimecompl(sprime,0.);
  std::complex<double> result = (sprime-(m1+m2)*(m1+m2))*PawianQFT::ChewMandelstamReid(sprimecompl, m1, m2);
    return result;
}

double CMIntegration2DReid::dsNorm(double sprime, double m1, double m2, double mpole, double fpole){
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


