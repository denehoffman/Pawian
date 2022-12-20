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

#include "gslUtils/CMIntegration2D.hh"
#include "Utils/Faddeeva.hh"
#include "qft++Extension/PawianUtils.hh"
#include "Utils/PawianConstants.hh"
#include "ErrLogger/ErrLogger.hh"

#include <complex>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

CMunstable_params CMIntegration2D::_CMunstable_params2(0.9033435, 0.59399735, 0.13957, 0.1349768, 0.13957);
std::complex<double> CMIntegration2D::_currentS2(1.,0.);
double CMIntegration2D::_currentsprimeprime(0.);
std::complex<double> CMIntegration2D::_currentOuterResult(0.,0.);
std::complex<double> CMIntegration2D::_currentInnerResult(0.,0.);
bool CMIntegration2D::_calcRealPartInner(true);

CMIntegration2D::CMIntegration2D(double mpole1, double fpole1, double m11, double m21,double mpole2, double fpole2, double m12, double m22) :
  CMIntegration(mpole1, fpole1, 0., m11, m21)
  {
  _integLowerBorder2=(_CMunstable_params2._m1+_CMunstable_params2._m2)*(_CMunstable_params2._m1+_CMunstable_params2._m2);
}

CMIntegration2D::~CMIntegration2D(){
}

void CMIntegration2D::setCMparams2(CMunstable_params& theParams){
  _CMunstable_params2._mPole=theParams._mPole;
  _CMunstable_params2._fPole=theParams._fPole;
  _CMunstable_params2._mu=theParams._mu;
  _CMunstable_params2._m1=theParams._m1;
  _CMunstable_params2._m2=theParams._m2;
  _integLowerBorder2=(_CMunstable_params2._m1+_CMunstable_params2._m2)*(_CMunstable_params2._m1+_CMunstable_params2._m2);
}

void CMIntegration2D::doFit2(gsl_function& F, double& result, double& resulterr, std::string fitName){
  const size_t n=1000;
  gsl_integration_workspace* wsp1=gsl_integration_workspace_alloc(n);
  int status = gsl_integration_qagiu(&F, _integLowerBorder2, _epsabs, _epsrel, _limit, wsp1, &result, &resulterr);
  if (status) {
    InfoMsg << fitName << " status: " << status << endmsg;
    if(status==GSL_EROUND){
      double epsabs=1.e-7;
      double epsrel=1.e-7;
      status = gsl_integration_qagiu(&F, _integLowerBorder2, epsabs, epsrel, _limit, wsp1, &result, &resulterr);
      InfoMsg << fitName << " status V1: " << status << endmsg;
      if(status==GSL_EROUND){
      double epsabs=1.e-6;
      double epsrel=1.e-6;
      status = gsl_integration_qagiu(&F, _integLowerBorder2, epsabs, epsrel, _limit, wsp1, &result, &resulterr);
      InfoMsg << fitName << " status V2: " << status << endmsg;      
      }
    }
  }
  gsl_integration_workspace_free (wsp1);
  //delete wsp1;
}

