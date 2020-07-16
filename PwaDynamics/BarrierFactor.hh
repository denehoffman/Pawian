
//************************************************************************//
//                                                                        //
//  Copyright 2014 Bertram Kopf (bertram@ep1.rub.de)                      //
//                 Julian Pychy (julian@ep1.rub.de)                       //
//                 - Ruhr-Universität Bochum                              //
//                                                                        //
//  This file is part of Pawian.                                          //
//                                                                        //
//  Pawian is free software: you can redistribute it and/or modify        //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or     //
//  (at your option) any later version.                                   //
//                                                                        //
//  Pawian is distributed in the hope that it will be useful,             //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of        //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         //
//  GNU General Public License for more details.                          //
//                                                                        //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.       //
//                                                                        //
//************************************************************************//

#pragma once

#include <complex>

class BarrierFactor
{

public:
   BarrierFactor(int l, std::complex<double> q0, double qR);
   std::complex<double> BlattWeisskopfRatio(std::complex<double> q);
   std::complex<double> BlattWeisskopf(std::complex<double> q);
   static std::complex<double> BlattWeisskopfRatio(int l, std::complex<double> q,
						   std::complex<double> q0, double qR);
  static std::complex<double> BlattWeisskopfTensorRatio(int l, std::complex<double> q,
							std::complex<double> q0, double qR);
   static std::complex<double> BlattWeisskopf(int l, std::complex<double> q, double qR);
   static std::complex<double> BlattWeisskopf(int l, std::complex<double> z);
   static std::complex<double> BlattWeisskopfTensor(int l, std::complex<double> q,
						    double qR);
   static const double qRDefault;

private:
   double _l;
   double _qR;
   std::complex<double> _q0;
   std::complex<double> _B0;
 };
