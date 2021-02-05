//************************************************************************//
//									  //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)			  //
//          	   - Ruhr-Universit??t Bochum 				  //
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

#include <algorithm>
#include "qft++Extension/AbsComplexFunction.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++Extension/PawianUtils.hh"

AbsComplexFunction::AbsComplexFunction()
{
}

AbsComplexFunction::~AbsComplexFunction(){
}

complex<double> AbsComplexFunction::eval(double s){
  return complex<double>(1.,0.);
}

complex<double> AbsComplexFunction::eval(complex<double> s){
  return complex<double>(1.,0.);
}

complex<double> AbsComplexFunction::interpolate(double current_s){
    std::vector<double>::iterator it;
    it = std::upper_bound(_s.begin(), _s.end(), current_s);
    int low = (it-_s.begin())-1;
    int up = (it-_s.begin());
    if(up == (int)_s.size()){
        low -= 1;
        up -= 1;
    }

    // get s, real and imag for lower and upper tabulated values
    double slow = _s.at(low); double sup = _s.at(up);
    double reallow = _realp.at(low); double realup = _realp.at(up);
    double imaglow = _imagp.at(low); double imagup = _imagp.at(up);

    // linear interpolation to actual s
    double m_real = (realup - reallow)/(sup - slow);
    double b_real = (realup + reallow-m_real*(sup + slow))/2.;
    double actual_real = m_real * current_s + b_real;

    double m_imag = (imagup - imaglow)/(sup - slow);
    double b_imag = (imagup + imaglow-m_imag*(sup + slow))/2.;
    double actual_imag = m_imag * current_s + b_imag;

    complex<double>amp = complex<double>(actual_real,actual_imag);
    return amp;
}


complex<double> AbsComplexFunction::interpolate(complex<double>current_s){
  return complex<double>(0,0);
}

