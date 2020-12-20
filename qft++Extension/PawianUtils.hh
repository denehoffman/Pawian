//************************************************************************//
//									  //
//  Copyright 2016 Bertram Kopf (bertram@ep1.rub.de)			  //
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

#pragma once

#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include <boost/math/special_functions/factorials.hpp>
#include <complex>

namespace PawianQFT{
/// Returns all valid @a LS combos for weak decays; non partity conservation
  vector<LS> GetValidLSWeak(const Spin &__j, const Spin &__s1, const Spin &__s2);
  complex<double> ChewMandelstamAsner(double s, double massDec1, double massDec2);
  double rhoiAsner(double s, double massDec1, double massDec2);
  complex<double> ChewMandelstamPennington(complex<double> s, double massDec1, double massDec2);
  complex<double> ChewMandelstamReid(complex<double> s, double massDec1, double massDec2);
  complex<double> ChewMandelstamReid_AngularMomentum(complex<double> s, double massDec1, double massDec2, int Spin);
  complex<double> ChewMandelstamDudek(complex<double> s, double massDec1, double massDec2);
  complex<double> ChewMandelstamDudekUnstableRhoPi(double s, double massDec1, double massDec2);
  complex<double> ChewMandelstamDudekUnstableRhoPi(complex<double> s, double massDec1, double massDec2);
  double evalHelpDudekUnstable(double mass, double mean, double sigma, double kl, double kh);
  complex<double> phaseSpaceFacReid(complex<double> mass, double massDec1, double massDec2);
  complex<double> phaseSpaceFacReid(double mass, double massDec1, double massDec2);
  complex<double> phaseSpaceFacDefault(double mass, double massDec1, double massDec2);
  complex<double> phaseSpaceFacDefault(complex<double> mass, double massDec1, double massDec2);
  complex<double> phaseSpaceFacAS(double mass, double massDec1, double massDec2);
  complex<double> phaseSpaceFacAS(complex<double> mass, double massDec1, double massDec2);
  double breakupMomQDefaultFromS(double s, double s1, double s2);
  template<typename MassType>
  complex<double> breakupMomQDefault(MassType mass, double massDec1, double massDec2);
  template<typename MassType>
  complex<double> breakupMomQReid(MassType mass, double massDec1, double massDec2);
  template<typename MassType>
  complex<double> breakupMomQDudek(MassType mass, double massDec1, double massDec2);
  template<typename MassType>
  complex<double> breakupMomQAS(MassType mass, double massDec1, double massDec2);
  template<typename MassType>
  complex<double> breakupMomQDefaultAS(MassType mass, double massDec1, double massDec2);
  complex<double> FlatteFkt(const Vector4<double> &__p4, std::pair<const double, const double>& decPair1, 
			    std::pair<const double, const double>& decPair2, double __mass0, double g1, double g2);
//
}
