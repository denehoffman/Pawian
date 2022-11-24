//************************************************************************//
//									  //
//  Copyright 2017 Markus Kuhlmann (mkuhlmann@ep1.rub.de)		  //
//  	      	   - Ruhr-Universität Bochum 				  //
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

#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "PwaDynamics/BBUnstableParFit.hh"
#include "PwaDynamics/BBUnstableParFitFcn.hh"
#include "qft++Extension/PawianUtils.hh"
#include <complex>


#include "ErrLogger/ErrLogger.hh"

BBUnstableParFit::BBUnstableParFit(double _mRes, double _wRes, double _m1, double _m2){
  mRes = _mRes;
  wRes = _wRes;
  m1 = _m1;
  m2 = _m2;
  // Display parameters for test distribution 

  cout << endl;
  InfoMsg <<"Set mRes as "<< mRes << endmsg;
  InfoMsg <<"Set wRes as "<< wRes << endmsg;
  InfoMsg <<"Set m1 as "<< m1 << endmsg;
  InfoMsg <<"Set m1 as "<< m2 << endmsg;
}

double BBUnstableParFit::calcNormD(const std::vector<double>& minPar){
  double result = 0.;
  std::complex<double> rootS(mRes,-wRes/2.0);
  std::complex<double> S = rootS*rootS;
  std::complex<double> intermediate;
  //intermediate = S-(minPar[0]*minPar[0])+(minPar[1]*minPar[1])*(S-((m1+m2)*(m1+m2)))*conj(PawianQFT::ChewMandelstamReid(mRes, m1, m2));
  intermediate = S-(minPar[0]*minPar[0])+(minPar[1]*minPar[1])*(S-((m1+m2)*(m1+m2)))*conj(PawianQFT::ChewMandelstamDudek(mRes, m1, m2));
  result = norm(intermediate);
  return result;
}

double BBUnstableParFit::calcNormDInvNeg(const std::vector<double>& minPar){
  double result = 0.;
  std::complex<double> rootS(mRes,-wRes/2.0);
  std::complex<double> S = rootS*rootS;
  std::complex<double> intermediate;
  //intermediate = S-(minPar[0]*minPar[0])+(minPar[1]*minPar[1])*(S-((m1+m2)*(m1+m2)))*conj(PawianQFT::ChewMandelstamReid(mRes, m1, m2));
  intermediate = S-(minPar[0]*minPar[0])+(minPar[1]*minPar[1])*(S-((m1+m2)*(m1+m2)))*conj(PawianQFT::ChewMandelstamDudek(mRes, m1, m2));
  result = -1.0/norm(intermediate);
  return result;
}

BBUnstableParFit::~BBUnstableParFit()
{
   
}

