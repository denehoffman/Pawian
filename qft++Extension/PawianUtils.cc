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
#include <boost/math/special_functions/factorials.hpp>
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

complex<double> PawianQFT::ChewMandelstamAsner(double s, double massDec1, double massDec2){
  double rho_i=rhoiAsner(s, massDec1, massDec2);
  complex<double> irhoComplex(-rho_i/PawianConstants::pi,0.);
  if(s<0 || s > (massDec1+massDec2)*(massDec1+massDec2)){
    double multTerm= fabs((1.+rho_i)/(1.-rho_i));
    irhoComplex*=log(multTerm);
    if(s > (massDec1+massDec2)*(massDec1+massDec2)) irhoComplex+=PawianConstants::i*rho_i; 
  }
  else{
    irhoComplex*=2.*atan(1/rho_i);
  }

  complex<double> result=irhoComplex;
  return result;  
}

double PawianQFT::rhoiAsner(double s, double massDec1, double massDec2){
  double resultSqr=fabs(1.-(massDec1+massDec2)*(massDec1+massDec2)/s);
  double result=sqrt(resultSqr);
  return result;  
}

complex<double> PawianQFT::ChewMandelstamPennington(complex<double> s, double massDec1, double massDec2){
  //Chew-Mandelstam parametrization
  //fulfils analyticity and unitarity
  // Basevant/Berger: ANL-HEP-PR-78-27 
  complex<double> sqrrho_a=complex<double>(1.,0.)-(massDec1+massDec2)*(massDec1+massDec2)/s;
  complex<double> rho_a = sqrt(sqrrho_a);
  complex<double> result=-rho_a/PawianConstants::pi*log((rho_a+1.)/(rho_a-1.));
  return result;
}

complex<double> PawianQFT::ChewMandelstamReid(complex<double> s, double massDec1, double massDec2){
  if(norm(s)<1.e-8) s=complex<double>(1.e-10, 1.e-10);
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

complex<double> PawianQFT::Power_Lambda(complex<double> s, double massDec1, double massDec2, double Power){
    complex<double> m1_p_m2_2 ((massDec1+massDec2)*(massDec1+massDec2), 0.);
    complex<double> m1_m_m2_2 ((massDec1-massDec2)*(massDec1-massDec2), 0.);

    complex<double>Lambda(0.,0.);

    Lambda += pow((m1_p_m2_2-s),Power)*pow((m1_m_m2_2-s),Power)/pow(s*s,Power);

    return Lambda;
}

complex<double> PawianQFT::ChewMandelstamReid_AngularMomentum(complex<double> s, double massDec1, double massDec2, int orbMom){

    // This Chew-Mandelstam Function is based on J. H. Reid, and N. N. Trofimenkoff, Journal of Mathematical Physics 25, 3540 (1984).
    // Here eq. (3.7) in case of equal masses and eq. (3.15) in case of unequal masses are implemented. 
    // In the s dependent log-Term the sign in the sqrt had to be changed to fulfill analyticity and unitarity which seemes to be a mistake in the paper.
    // All output of this function has been checked carefully and fullfils analyticity and unitarity.
    // Meike
    
    if(norm(s)<1.e-8) s=complex<double>(1.e-10, 1.e-10);

    if(abs(massDec1 - massDec2)<1.e-8){

        complex<double> k2 = (s-pow(massDec1+massDec2,2.))*(s-pow(massDec1-massDec2,2.))/(4.*s);
        complex<double> massDec(massDec1,0.);
        complex<double> lambda = (4.*k2)/s;
        complex<double> sum (0.,0.);

        for(int j = 0; j<=orbMom; j++){
            sum+=pow(lambda, orbMom-j)*1./(2.*j+1.);
        }

        complex<double> CM = 2./PawianConstants::pi*(-pow(lambda, ((double)orbMom+0.5))*log((sqrt(-s+4.*massDec*massDec)+sqrt(-s))/(2.*massDec))+sum);

        complex<double> result = CM;
        
        return result;
    }

    else{
        
        if(norm(s)<1.e-8) s=complex<double>(1.e-10, 1.e-10);

        complex<double> CM (0.,0.);

        if(norm(pow(s,orbMom))<1.e-10) return CM;

        complex<double> m1_p_m2_2((massDec1+massDec2)*(massDec1+massDec2), 0.);
        complex<double> m1_m_m2_2((massDec1-massDec2)*(massDec1-massDec2), 0.);

        complex<double> m1_p_m2_2_m_s = (m1_p_m2_2-s);
        complex<double> m1_m_m2_2_m_s = (m1_m_m2_2-s);

        complex<double> k2 = (s-m1_p_m2_2)*(s-m1_m_m2_2)/(4.*s);
        complex<double> lambda = (4.*k2)/s;
        
        complex<double> mu = pow((m1_p_m2_2-m1_m_m2_2),2.)/(16.*m1_p_m2_2*m1_m_m2_2);
        complex<double> nu = (m1_p_m2_2+m1_m_m2_2)/(2.*sqrt(m1_p_m2_2*m1_m_m2_2));
        //complex<double> omega = (massDec1*massDec1+massDec2*massDec2)/(massDec1*massDec1-massDec2*massDec2)-(massDec1*massDec1-massDec2*massDec2)/s; 

        complex<double> sum1 (0.,0.); 
        complex<double> sum2 (0.,0.);
        complex<double> sum3 (0.,0.);

        double Omega=0.;

        for(int j = 0; j<=orbMom; j++){
            
            Omega=pow(-1.,j)*pow(2.,j)*boost::math::double_factorial<double>(2*j-1)/boost::math::factorial<double>(j);
            sum1 += PawianQFT::Power_Lambda(s, massDec1, massDec2, orbMom-j)/(2.*j+1.);
            sum2 += Omega*PawianQFT::Power_Lambda(s, massDec1, massDec2, orbMom-j)*pow(mu,j);
        }

        for(int p = 1; p<=orbMom; p++){

            complex<double> sumA (0.,0.);
            complex<double> sumB (0.,0.);

            for(int q = 0; q<=orbMom-p; q++){
                Omega=pow(-1.,q)*pow(2.,q)*boost::math::double_factorial<double>(2*q-1)/boost::math::factorial<double>(q);
                sumA+=Omega*Power_Lambda(s, massDec1, massDec2, orbMom-q-p)*pow(mu,q);
            }

            for(int r = 0; r<=p-1; r++){
                Omega=pow(-1.,r)*pow(2.,r)*boost::math::double_factorial<double>(2*r-1)/boost::math::factorial<double>(r);
                sumB+=Omega*pow(mu,r);
            }

            sum3+=sumA*sumB*(1./p);

        }

        CM += -2./PawianConstants::pi*(-sqrt(m1_p_m2_2_m_s)*sqrt(m1_m_m2_2_m_s)/s*pow(lambda,orbMom)*log((sqrt(m1_p_m2_2-s)+sqrt(m1_m_m2_2-s))/(2.*sqrt(massDec1*massDec2)))-0.5*sum1 - 0.5*(massDec1*massDec1+massDec2*massDec2)/(massDec1*massDec1-massDec2*massDec2)*log(massDec1/massDec2)*sum2 + 0.5*(massDec1*massDec1-massDec2*massDec2)/s*log(massDec1/massDec2)*sum2 - nu*nu*sum3/4. + sqrt(m1_p_m2_2*m1_m_m2_2)/s*nu*sum3/4.);

        return CM;
    }

}

double PawianQFT::Offset_Dudek_Reid(double massDec1, double massDec2){

    double result = 0.;

    if(abs(massDec1-massDec2)<1e-6){
        result = (1./PawianConstants::pi*(pow(massDec1,2.)-pow(massDec2,2.) + 2.* massDec1*massDec2*log(massDec1/massDec2)))/(massDec1*massDec1-massDec2*massDec2);
    } else {
        result = 2./PawianConstants::pi;
    }

    return result;
}

complex<double> PawianQFT::ChewMandelstamDudek_AngularMomentum(complex<double> s, double massDec1, double massDec2, int orbMom){

    complex<double> CM_ReidMom = PawianQFT::ChewMandelstamReid_AngularMomentum(s, massDec1, massDec2, orbMom);
    double Offset = PawianQFT::Offset_Dudek_Reid(massDec1, massDec2);

    complex<double> result(Offset-CM_ReidMom.real(),-CM_ReidMom.imag());
    return result;
}

complex<double> PawianQFT::ChewMandelstamDudek(complex<double> s, double massDec1, double massDec2){
  if(norm(s)<1.e-8) s=complex<double>(1.e-10, 1.e-10);
  complex<double> m1_2_m_m2_2(massDec1*massDec1-massDec2*massDec2, 0.);
  complex<double> m1_2_p_m2_2(massDec1*massDec1+massDec2*massDec2, 0.); 

  complex<double> m1_p_m2_2((massDec1+massDec2)*(massDec1+massDec2), 0.);
  complex<double> m1_m_m2_2((massDec1-massDec2)*(massDec1-massDec2), 0.);

  if(abs(s - m1_p_m2_2)<1.e-10) return complex<double>(0.,0.);

  complex<double> rho(0.,0.);
  rho += sqrt(complex<double>(1,0)-(m1_p_m2_2/s))*sqrt(complex<double>(1,0)-(m1_m_m2_2/s));

  complex<double> xi(1.,0.);
  xi -= (m1_p_m2_2/s);


  complex<double> term1(1., 0.);
  term1 *= rho/PawianConstants::pi * log((xi+rho)/(xi-rho));

  complex<double> term2(1., 0.);
  term2 *= xi/PawianConstants::pi * (massDec2-massDec1)/(massDec1+massDec2) * log(massDec2/massDec1);


  complex<double> result = term1 - term2;
  return result;
}

complex<double> PawianQFT::ChewMandelstamDudekUnstableRhoPi(double s, double massDec1, double massDec2){
  double sqrts=sqrt(s);
//  double parImag[17]={-56.9508, 257.776, -499.825, 540.399, 0.510786, 0.417175, 0.091638, -1.28979, 7.27124, -21.048, -352.435, 137.719, -26.7861, 33.0017, -26.5003, 8.37052, 0.279228};
//  double parReal[12]={0.913953, -1.63788, 1.2278, -0.602585, 0.945479, 0.00101412, 0.377936, -0.361194, 0.332049, -0.0618304, 0.00291903, 0.0124467};

//     double resultReal=(parReal[0]+parReal[1]*sqrts+parReal[2]*sqrts*sqrts+parReal[3]*pow(sqrts,3))*evalHelpDudekUnstable(sqrt(sqrts), parReal[4], parReal[5], parReal[10], parReal[11])+parReal[6]+parReal[7]*sqrts+parReal[8]*sqrts*sqrts+parReal[9]*pow(sqrts,3.);

//  double resultReal=(parReal[0]+parReal[1]*sqrts+parReal[2]*sqrts*sqrts+parReal[3]*pow(sqrts,3))*evalHelpDudekUnstable(sqrt(sqrts), parReal[4], parReal[5], parReal[10], parReal[11])+parReal[6]+parReal[7]*sqrts+parReal[8]*sqrts*sqrts+parReal[9]*pow(sqrts,3.);

//     double resultImag=(parImag[0]+parImag[1]*sqrts+parImag[2]*sqrts*sqrts+parImag[3]*sqrts*sqrts*sqrts+parImag[10]*sqrts*sqrts*sqrts*sqrts+parImag[11]*pow(sqrts,5.)+parImag[12]*pow(sqrts,6.))*phaseSpaceFacDefault(sqrts, parImag[4], parImag[5]).imag()+parImag[6]+parImag[7]*sqrts+parImag[8]*sqrts*sqrts+parImag[9]*pow(sqrts,3.)+parImag[13]*pow(sqrts,4.)+parImag[14]*pow(sqrts,5.)+parImag[15]*pow(sqrts,6.)+parImag[16]*pow(sqrts,8.);

  double parImag[18]={ -56.62759611094936, 260.0512366895287, -520.6462126261961, 597.0275352864451, 0.5651508581192156, 0.345639612989795, -2.137663230250565, -30.11215274793011, 56.77076445863649, -98.83065949201345, -426.9460849355088, 187.1610397574648, -40.14298412560645, 110.1522287917516, -66.53138741774561, 13.87591536254525, 2.41492594582418, 14.58588929359466};
  
  double parReal[14]={1.993577947915189, -2.904436961229936, -0.1390783363995473, 1.96131811547056, 0.9207725124387273, 0.1220254365759813, -2.079685061379917, 5.1232003549247, -3.766588111727992, 1.563541895855885, -0.04030295032523248, 0.5840845557110333, -1.349512674240738, -0.2487791849975481};

  double resultReal=(parReal[0]+parReal[1]*sqrts+parReal[2]*sqrts*sqrts+parReal[3]*pow(sqrts,3)+parReal[12]*pow(sqrts,4))*evalHelpDudekUnstable(sqrt(sqrts), parReal[4], parReal[5], parReal[10], parReal[11])+parReal[6]+parReal[7]*sqrts+parReal[8]*sqrts*sqrts+parReal[9]*pow(sqrts,3.)+parReal[13]*pow(sqrts,4);


double resultImag = (parImag[0]+parImag[1]*sqrts+parImag[2]*sqrts*sqrts+parImag[3]*sqrts*sqrts*sqrts+parImag[10]*sqrts*sqrts*sqrts*sqrts+parImag[11]*pow(sqrts,5.)+parImag[12]*pow(sqrts,6.))*phaseSpaceFacDefault(sqrts, parImag[4], parImag[5]).real()+parImag[6]+parImag[7]*sqrts+parImag[8]*sqrts*sqrts+parImag[9]*pow(sqrts,3.)+parImag[13]*pow(sqrts,4.)+parImag[14]*pow(sqrts,5.)+parImag[15]*pow(sqrts,6.)+parImag[16]*pow(sqrts,7.)+parImag[17]*sqrt(sqrts);

     complex<double> result(resultReal, -resultImag);
     return result;
}

double PawianQFT::evalHelpDudekUnstable(double mass, double mean, double sigma, double kl, double kh){
  double val=0.;
  double arg=(mass-mean)/sigma;

  if(arg<=-kl){
    val=exp(0.5*kl*kl+kl*arg);
  }
  
  if(arg>-kl&&arg<=kh){
    val=exp(-0.5*arg*arg);
  }
  
  if(arg>kh){
    val=exp(0.5*kh*kh-kh*arg);
  }
  
  return val;
}

complex<double> PawianQFT::ChewMandelstamDudekUnstableRhoPi(complex<double> s, double massDec1, double massDec2){
  complex<double> result(0.,0.);
  return result; //dummy
}

complex<double> PawianQFT::phaseSpaceFacReid(complex<double> mass, double massDec1, double massDec2){
  complex<double> s=mass*mass;
  complex<double> result=PawianQFT::ChewMandelstamReid(s, massDec1, massDec2).imag();
  return result; 
}

complex<double> PawianQFT::phaseSpaceFacReid(double mass, double massDec1, double massDec2){
  complex<double> massCompl(mass, 1.e-10); // for real s: expansion to s=0 from 1st quadrant
  return PawianQFT::phaseSpaceFacReid(massCompl, massDec1, massDec2);
} 

complex<double> PawianQFT::phaseSpaceFacDefault(double mass, double massDec1, double massDec2){

  complex<double> result(0.,0.);

  if(fabs(mass) < 1e-8) {
     std::cout << "mass " << mass 
	       << " very close to 0; not possible to calculate phasespace factor: \n set mass to 1.*e-10 !!!" 
	       << std::endl;
     mass=1.e-10;
  }

  double termPlus=(massDec1+massDec2)/mass;
  double termMinus=(massDec1-massDec2)/mass;
  double tmpVal=(1.-termPlus*termPlus) * (1.-termMinus*termMinus);
  if(tmpVal>=0.) result = std::complex<double>(std::sqrt(tmpVal), 0.);
  else result = std::complex<double>(0., std::sqrt(-tmpVal));   
  return result;
}

complex<double> PawianQFT::phaseSpaceFacDefault(complex<double> mass, double massDec1, double massDec2){

  complex<double> result(0.,0.);

  if(norm(mass) < 1e-8) {
     std::cout << "mass " << mass 
	       << " very close to 0; not possible to calculate phasespace factor: \n set real part of the mass to 1.*e-10 !!!" 
	       << std::endl;
     mass=complex<double>(1.e-10,0.);
  }

  complex<double> termPlus=(massDec1+massDec2)/mass;
  complex<double> termMinus=(massDec1-massDec2)/mass;
  complex<double> tmpVal=(1.-termPlus*termPlus) * (1.-termMinus*termMinus);
  result = std::sqrt(tmpVal);
  return result;
}

complex<double> PawianQFT::phaseSpaceFacAS(double mass, double massDec1, double massDec2){

  complex<double> result(0.,0.);

  if(fabs(mass) < 1e-8) {
     std::cout << "mass " << mass 
	       << " very close to 0; not possible to calculate phasespace factor: \n set mass to 1.*e-10 !!!" 
	       << std::endl;
     mass=1.e-10;
  }

  double termPlus=(massDec1+massDec2)/mass;
  double tmpVal=(1.-termPlus*termPlus);
  if(tmpVal>=0.) result = std::complex<double>(std::sqrt(tmpVal), 0.);
  else result = std::complex<double>(0., std::sqrt(-tmpVal));   
  return result;
}

complex<double> PawianQFT::phaseSpaceFacAS(complex<double> mass, double massDec1, double massDec2){

  complex<double> result(0.,0.);

  if(norm(mass) < 1e-8) {
     std::cout << "mass " << mass 
	       << " very close to 0; not possible to calculate phasespace factor: \n set real part of the mass to 1.*e-10 !!!" 
	       << std::endl;
     mass=complex<double>(1.e-10,0.);
  }

  complex<double> termPlus=(massDec1+massDec2)/mass;
  complex<double> tmpVal=(1.-termPlus*termPlus);
  result = std::sqrt(tmpVal);
  return result;
}

double PawianQFT::breakupMomQDefaultFromS(double s, double s1, double s2){
  if ((s*s + s1*s1 + s2*s2 - 2. * ( s*s1 + s*s2 + s1*s2)) < 0.){
    std::cout << "(s*s + s1*s1 + s2*s2 - 2. * ( s*s1 + s*s2 + s1*s2)) = "
	      << (s*s + s1*s1 + s2*s2 - 2. * ( s*s1 + s*s2 + s1*s2))
              << " < 0  -> caluclation is not allowed!!!!"
	      << std::endl;
    exit(1);
  }
  double result = std::sqrt(s*s + s1*s1 + s2*s2 - 2. * ( s*s1 + s*s2 + s1*s2))/std::sqrt(4.*s);
  return result;
}


template<typename MassType>
complex<double> PawianQFT::breakupMomQDefault(MassType mass, double massDec1, double massDec2){
  complex<double> result=PawianQFT::phaseSpaceFacDefault(mass, massDec1, massDec2)*mass/2.;
  return result;  
}

template<typename MassType>
complex<double> PawianQFT::breakupMomQReid(MassType mass, double massDec1, double massDec2){
  double rho = PawianQFT::ChewMandelstamReid(mass*mass, massDec1, massDec2).imag();
  //  if(rho < 0) rho*=-1.; 
  complex<double> result=rho*mass/2.;
  return result;  
}

template<typename MassType>
complex<double> PawianQFT::breakupMomQDudek(MassType mass, double massDec1, double massDec2){
  double rho = -PawianQFT::ChewMandelstamDudek(mass*mass, massDec1, massDec2).imag();
  complex<double> result=rho*mass/2.;
  return result;  
}

template<typename MassType>
complex<double> PawianQFT::breakupMomQAS(MassType mass, double massDec1, double massDec2){
  complex<double> result=PawianQFT::phaseSpaceFacAS(mass, massDec1, massDec2)*mass/2.;
  return result;  
}

template<typename MassType>
complex<double> PawianQFT::breakupMomQDefaultAS(MassType mass, double massDec1, double massDec2){
  complex<double> result=PawianQFT::breakupMomQDefault(mass, massDec1, massDec2);
  complex<double> qAS=PawianQFT::breakupMomQAS(mass, massDec1, massDec2);
  if( fabs(qAS.real()) < 1.e-10 ) result=qAS;
  return result; 
}

template complex<double> PawianQFT::breakupMomQDefault(double, double, double);
template complex<double> PawianQFT::breakupMomQDefault(complex<double>, double, double);
template complex<double> PawianQFT::breakupMomQReid(double, double, double);
template complex<double> PawianQFT::breakupMomQReid(complex<double>, double, double);
template complex<double> PawianQFT::breakupMomQDudek(double, double, double);
template complex<double> PawianQFT::breakupMomQDudek(complex<double>, double, double);
template complex<double> PawianQFT::breakupMomQAS(double, double, double);
template complex<double> PawianQFT::breakupMomQAS(complex<double>, double, double);
template complex<double> PawianQFT::breakupMomQDefaultAS(double, double, double);
template complex<double> PawianQFT::breakupMomQDefaultAS(complex<double>, double, double);

complex<double> PawianQFT::FlatteFkt(const Vector4<double> &__p4, std::pair<const double, const double>& decPair1, std::pair<const double
		       , const double>& decPair2, double __mass0, double g1, double g2){

  complex<double> i(0.,1.);
  const double m1a=decPair1.first;
  const double m1b=decPair1.second;
  const double m2a=decPair2.first;
  const double m2b=decPair2.second;

  double mAB=__p4.Mass();

  //calculate gammas with phase-space factors 
  complex<double> gamma11=g1*PawianQFT::breakupMomQDefault(mAB, m1a, m1b);
  complex<double> gamma22=g2*PawianQFT::breakupMomQDefault(mAB, m2a, m2b);

  complex<double> gammaLow(0.,0.);
  if( (m1a+m1b) < (m2a+m2b) ) gammaLow=gamma11;
  else gammaLow=gamma22;

  complex<double>  result=__mass0*sqrt(gammaLow*gamma11)/( __mass0*__mass0 - mAB*mAB - i * __mass0 * (gamma11+gamma22) );

  return result;

}

