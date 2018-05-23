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

#pragma once 

#include <iostream>
#include <vector>
#include "PwaDynamics/PhaseSpaceIsobar.hh"
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

using namespace std;

//_____________________________________________________________________________
//_____________________________________________________________________________
struct LUTLine {
  complex<double> requested;
  complex<double> cm;
} ;

class PhaseSpaceIsobarLUT: public PhaseSpaceIsobar {

public:
  
  /// Constructor
  PhaseSpaceIsobarLUT(double mass1, double mass2, string type); 

  /// Destructor
  virtual ~PhaseSpaceIsobarLUT();

  // operators:


  // functions:

  virtual complex<double> factor(const double mass);
  virtual complex<double> breakUpMom(const double mass);
  virtual complex<double> factor(const complex<double> mass);
  virtual complex<double> breakUpMom(const complex<double> mass);
  virtual complex<double> ChewM(const double mass);
  virtual complex<double> ChewM(const complex<double> mass);
  virtual void cacheFactors(const double mass);
  complex<double> sLow(){
    return m_sLow;
  };
  complex<double> sHigh(){
    return m_sHigh;
  };

protected:
  void CorrectCMForChosenSign(complex<double>& breakUpMom, complex<double>& toChange);

private:
  LUTLine readLine(int);
  complex<double> interpolateCM(LUTLine, LUTLine, LUTLine, LUTLine, complex<double>);
  complex<double> interpolateCM(LUTLine, LUTLine);
  std::ifstream m_lutfile;
  std::string m_lutfilepath;
  complex<double> getFactor(complex<double>);
  int m_dSize;
  int m_nRe;
  int m_nIm;
  complex<double> m_sLow;
  complex<double> m_sHigh;
  void loadParams();
  bool m_sortedByReS;
};
//

