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
#include "ErrLogger/ErrLogger.hh"

#include <complex>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

CMunstable_params CMIntegration::_CMunstable_params(0.9033435, 0.59399735, 0.13957, 0.1349768, 0.13957);
bool CMIntegration::_calcRealPart=true;
std::complex<double> CMIntegration::_currentS(1.,0.);

CMIntegration::CMIntegration(double mpole, double fpole, double mu, double m1, double m2):
  _epsabs(1e-8)
  ,_epsrel(1e-8)
  ,_limit(1000){
  //_CMunstable_params._mPole=mpole;
  gsl_set_error_handler_off();
  _integLowerBorder=(_CMunstable_params._m1+_CMunstable_params._m2)*(_CMunstable_params._m1+_CMunstable_params._m2);
}

CMIntegration::~CMIntegration(){
}

void CMIntegration::setCMparams(CMunstable_params& theParams){
  _CMunstable_params._mPole=theParams._mPole;
  _CMunstable_params._fPole=theParams._fPole;
  _CMunstable_params._mu=theParams._mu;
  _CMunstable_params._m1=theParams._m1;
  _CMunstable_params._m2=theParams._m2;
}

void CMIntegration::doFit(gsl_function& F, double& result, double& resulterr, std::string fitName){
  const size_t n=1000;
  gsl_integration_workspace* wsp1=gsl_integration_workspace_alloc(n);
  int status = gsl_integration_qagiu(&F, _integLowerBorder, _epsabs, _epsrel, _limit, wsp1, &result, &resulterr);
  if (status) {
    InfoMsg << fitName << " status: " << status << endmsg;
    if(status==GSL_EROUND){
      double epsabs=1.e-7;
      double epsrel=1.e-7;
      status = gsl_integration_qagiu(&F, _integLowerBorder, epsabs, epsrel, _limit, wsp1, &result, &resulterr);
      InfoMsg << fitName << " status V1: " << status << endmsg;
      if(status==GSL_EROUND){
      double epsabs=1.e-6;
      double epsrel=1.e-6;
      status = gsl_integration_qagiu(&F, _integLowerBorder, epsabs, epsrel, _limit, wsp1, &result, &resulterr);
      InfoMsg << fitName << " status V2: " << status << endmsg;      
      }
    }
  }
}

