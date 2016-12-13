//************************************************************************//
//                                                                        //
//  Copyright 2016 Bertram Kopf (bertram@ep1.rub.de)                      //
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

// DataUtils class definition file. -*- C++ -*-
// Copyright 2016 Bertram Kopf

#include "qft++Extension/PawianUtils.hh"
#include "Utils/PawianConstants.hh"
#include <math.h>

vector<LS> PawianQFT::GetValidLSWeak(const Spin &__j, const Spin &__s1, const Spin &__s2){
  vector<LS> valid_ls;
  LS ls;

  for(Spin S = abs(__s1 - __s2); S <= (__s1 + __s2); ++S){
    for(int L = (int)abs(__j - S); L <= (int)(__j + S); ++L){
        ls.L = L;
        ls.S = S;
        valid_ls.push_back(ls);
    }
  }
  return valid_ls;
}

complex<double> PawianQFT::phaseSpaceFacAsner(double s, double massDec1, double massDec2){
  complex<double> i(0.,1.);

  double rho_i=rhoiAsner(s, massDec1, massDec2);
  complex<double> irhoComplex(-rho_i/PawianConstants::pi,0.);
  if(s<0 || s > (massDec1+massDec2)*(massDec1+massDec2)){
    double multTerm= fabs((1.+rho_i)/(1.-rho_i));
    irhoComplex*=log(multTerm);
    if(s > (massDec1+massDec2)*(massDec1+massDec2)) irhoComplex+=i*rho_i; 
  }
  else{
    irhoComplex*=2.*atan(1/rho_i);
  }

  //  complex<double> result=irhoComplex/i;
  complex<double> result=irhoComplex;
  return result;  
}

double PawianQFT::rhoiAsner(double s, double massDec1, double massDec2){
  double resultSqr=fabs(1.-(massDec1+massDec2)*(massDec1+massDec2)/s);
  double result=sqrt(resultSqr);
  return result;  
}

complex<double> PawianQFT::phaseSpaceFacPennington(complex<double> s, double massDec1, double massDec2){
  //Chew-Mandelstam parametrization
  //fulfils analyticity and unitarity
  // Basevant/Berger: ANL-HEP-PR-78-27 
 complex<double> i(0.,1.);
  complex<double> sqrrho_a=complex<double>(1.,0.)-(massDec1+massDec2)*(massDec1+massDec2)/s;
  complex<double> rho_a = sqrt(sqrrho_a);
  complex<double> result=-rho_a/PawianConstants::pi*log((rho_a+1.)/(rho_a-1.));
  return result;
}

complex<double> PawianQFT::ChewMandelstamReid(complex<double> s, double massDec1, double massDec2){
  if(norm(s)<1.e-8) s=complex<double>(1.e-10,0.);
  complex<double> m1_2_m_m2_2(massDec1*massDec1-massDec2*massDec2, 0.);
  complex<double> m1_2_p_m2_2(massDec1*massDec1+massDec2*massDec2, 0.); 

  complex<double> m1_p_m2_2((massDec1+massDec2)*(massDec1+massDec2), 0.);
  complex<double> m1_m_m2_2((massDec1-massDec2)*(massDec1-massDec2), 0.);

  complex<double> m1_p_m2_2_m_s=m1_p_m2_2-s;
  complex<double> m1_m_m2_2_m_s=m1_m_m2_2-s;

  complex<double> term1(1., 0.);
  term1*=sqrt(m1_p_m2_2_m_s)*sqrt(m1_m_m2_2_m_s)*log((sqrt(m1_p_m2_2_m_s)+sqrt( m1_m_m2_2_m_s))/(2.*sqrt(massDec1*massDec2)))/s;


  complex<double> term2(0.5, 0.);
  term2*=m1_2_m_m2_2*log(massDec1/massDec2)/s;


  complex<double> term34(0.5, 0.);
  if( norm(m1_2_m_m2_2) > 1.e-15) term34+=log(massDec1/massDec2)*m1_2_p_m2_2/(2.*m1_2_m_m2_2);
  else term34+=complex<double>(0.5, 0.);  

  complex<double> result=-2./PawianConstants::pi*(term2-term1-term34);
  return result;
}

complex<double> PawianQFT::phaseSpaceFacReid(complex<double> mass, double massDec1, double massDec2){
  complex<double> i(0.,1);
  complex<double> s=mass*mass;
  complex<double> result=-PawianQFT::ChewMandelstamReid(s, massDec1, massDec2)/i;
  return result; 
}

complex<double> PawianQFT::phaseSpaceFacReid(double mass, double massDec1, double massDec2){
  complex<double> massCompl(mass, 0.);
  return PawianQFT::phaseSpaceFacReid(massCompl, massDec1, massDec2); 
}

complex<double> PawianQFT::phaseSpaceFacPenningtonsqrts(complex<double> sqrts, double massDec1, double massDec2){
  complex<double> s=sqrts*sqrts;
  complex<double> result=PawianQFT::phaseSpaceFacPennington(s, massDec1, massDec2);
  return result;
}


