#include <algorithm>   
#include <cmath>

#include "qft++Extension/SphericalHarmonic.hh"
#include "Utils/MathUtils.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"

namespace SphericalHarmonic{
  double LegendrePolynomial(unsigned int l, int m, double cosT){
    int LplusM=l+m;
    double prefac1= MathUtils::factorial(LplusM);

    int LminusM=l-m;
    double prefac2= MathUtils::factorial(LminusM);

    double result=prefac1/prefac2*Wigner_d((int) l, m, 0, acos(cosT));    
    return result;
  }

  std::complex<double> SphericalHarmonicFunc(int l, int m, double theta, double phi){
    // double prefac1=(2.*l+1.)/(4*M_PI);
    // int LminusM=l-m;
    // double prefac2= MathUtils::factorial(LminusM);
    // int LplusM=l+m;
    // double prefac3= MathUtils::factorial(LplusM);
    // double prefacAll=sqrt(prefac1*prefac2/prefac3);
    // //    if(m>0) prefacAll*=pow(-1.,m);
    
    // std::complex<double> exp_imphi(cos(m*phi), sin(m*phi));
    // std::complex<double> result=prefacAll*LegendrePolynomial(l, m, cos(theta))*exp_imphi;

    double prefac=sqrt((2.*l+1.)/(4*M_PI));
    std::complex<double> result=prefac*conj(Wigner_D(phi, theta, 0, l, m, 0));
    if(result.imag()!=0.){
      std::cout << "result: " << result << std::endl;
    }    
    return result;
  }
}
