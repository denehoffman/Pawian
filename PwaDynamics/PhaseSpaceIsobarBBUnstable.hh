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
#include <julia.h>
#include <map>
#include <string>
#include <sstream>

using namespace std;

//_____________________________________________________________________________
//_____________________________________________________________________________


class PhaseSpaceIsobarBBUnstable: public PhaseSpaceIsobar {

public:
  
  /// Constructor
  PhaseSpaceIsobarBBUnstable(double mass1, double mass2, string type); 

  /// Destructor
  virtual ~PhaseSpaceIsobarBBUnstable();

  // operators:


  // functions:

  virtual complex<double> factor(const double mass);
  virtual complex<double> breakUpMom(const double mass);
  virtual complex<double> factor(const complex<double> mass);
  virtual complex<double> breakUpMom(const complex<double> mass);
  virtual complex<double> ChewM(const double mass);
  virtual complex<double> ChewM(const complex<double> mass);
   virtual void cacheFactors(const double mass);

protected:

private:
  jl_module_t* m_mod;
  jl_function_t* m_func;
  jl_value_t* m_srealHO;
  jl_value_t* m_simagHO;
  jl_value_t* m_m1_1HO;
  jl_value_t* m_m1_2HO;
  jl_value_t* m_m2_1HO;
  jl_value_t* m_m2_2HO;
  jl_value_t* m_mR1HO;
  jl_value_t* m_mR2HO;
  jl_value_t* m_f1HO;
  jl_value_t* m_f2HO;
  jl_value_t* m_epsilonHO;

  std::string m_pathToModule;
  std::map<std::string, double*> m_knownCombinations;

  bool loadModule(std::string _moduleDefinitionPath);
  double* computeFactor(double* _inS);
  double* computeFactor(double _inSreal);
  complex<double> computeFactor(complex<double> inS);
  
  double* computeFactor();
  bool setArgs(double* _inArgs);
  bool setS(double* _inS);
  void fillMap();
};
//_____________________________________________________________________________


