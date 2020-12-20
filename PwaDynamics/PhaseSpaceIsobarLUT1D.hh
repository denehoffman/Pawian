//************************************************************************//
//									                                      //
//  Copyright 2020 Meike Kuessner (mkuessner@ep1.rub.de)		          //
//  	      	   - Ruhr-Universität Bochum 				              //
//									                                      //
//  This file is part of Pawian.			                    		  //
//									                                      //
//  Pawian is free software: you can redistribute it and/or modify	      //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	              //
//									                                      //
//  Pawian is distributed in the hope that it will be useful,		      //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	      //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	      //
//  GNU General Public License for more details.	      		          //
//									                                      //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	      //
//                                  									  //
//************************************************************************//

#pragma once 

#include "PwaDynamics/PhaseSpaceIsobar.hh"

#include <string>
#include <fstream>
#include <complex>
using std::complex;

struct LUTLine1D {
  double requested;
  complex<double> cm;
} ;

class PhaseSpaceIsobarLUT1D: public PhaseSpaceIsobar {

public:
  /// Constructor
  PhaseSpaceIsobarLUT1D(double mass1, double mass2, string type); 

  /// Destructor
  virtual ~PhaseSpaceIsobarLUT1D();

  virtual complex<double> factor(const double mass);
  virtual complex<double> factor(const complex<double> mass);
  virtual complex<double> breakUpMom(const double mass);
  virtual complex<double> breakUpMom(const complex<double> mass);
  virtual complex<double> ChewM(const double mass, int orbMom=0);
  virtual complex<double> ChewM(const complex<double> mass, int orbMom=0);
  virtual void cacheFactors(const double mass);

  double sLow() {
    return m_sLow;
  };
  double sHigh() {
    return m_sHigh;
  };

protected:
  void CorrectCMForChosenSign(complex<double>& breakUpMom, complex<double>& toChange);

private:
  LUTLine1D readLine(int);
  complex<double> interpolateCM(LUTLine1D, LUTLine1D, double);
  std::ifstream m_lutfile;
  std::string m_lutfilepath;
  complex<double> getFactor(double);
  int m_dSize;
  int m_nRe;
  int m_nIm;
  double m_sLow;
  double m_sHigh;
  void loadParams();
  bool m_sortedByReS;
};

