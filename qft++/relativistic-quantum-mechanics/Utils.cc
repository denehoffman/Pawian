#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++/relativistic-quantum-mechanics/BlattWeisskopf.hh"

/* Copyright 2008 Mike Williams (mwill@jlab.org)
 *
 * This file is part of qft++.
 *
 * qft++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * qft++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with qft++.  If not, see <http://www.gnu.org/licenses/>.
 */
//_____________________________________________________________________________
/** @file Utils.C
 *  @brief Utility function source file
 */
//_____________________________________________________________________________
/// Utility function used by Clebsch 
inline double dfact(double __x){
  if((__x < 0.00001) && (__x >= 0.0)) return 1.;
  if(__x < 0) return 0.;
  return __x*dfact(__x - 1.);
}
//_____________________________________________________________________________
/// Returns i!
inline int factorial(int __i) {
  if (__i<0){
    cerr << endl;
    cerr << "factoral value " << __i << " must be >=0 !!! "  
	 << endl;

     exit(1);
  }
  int f = 1;
  if((__i == 0) || (__i == 1)) f = 1;  
  else{
    while(__i > 0){
      f = f*__i;
      __i--;
    }
  }
  return f;
} 
//_____________________________________________________________________________
// define a few macros
#define MAX(x,y) (x>y ? x : y)
#define MIN(x,y) (x<y ? x : y)
//_____________________________________________________________________________

vector<LS> GetValidLS(const Spin &__j,int __parity,const Spin &__s1,int __p1,
		      const Spin &__s2,int __p2){

  vector<LS> valid_ls;
  LS ls;

  for(Spin S = abs(__s1 - __s2); S <= (__s1 + __s2); S++){
    for(int L = (int)abs(__j - S); L <= (int)(__j + S); L++){
      // check the parity
      if(abs(__p1*__p2*pow(-1.,L) - __parity) < 1.e-5) {
	ls.L = L;
	ls.S = S;
	valid_ls.push_back(ls);
      }
    }
  }
  return valid_ls;
}
//_____________________________________________________________________________

double Gamma(double __z){
  
  float x = abs(__z);
  // Coefficients for the series expansion
  double c[7] = { 2.5066282746310005, 76.18009172947146, -86.50532032941677
	       ,24.01409824083091, -1.231739572450155, 0.1208650973866179e-2
               ,-0.5395239384953e-5};

   double y = x;
   double tmp = x + 5.5;
   tmp = (x+0.5)*log(tmp)-tmp;
   double ser = 1.000000000190015;
   for(int i = 1; i < 7; i++) {
      y += 1;
      ser += c[i]/y;
   }
   double lnGamma = tmp+log(c[0]*ser/x); // natural log of Gamma(z)
   if(__z > 0) return exp(lnGamma);   
   else if(__z == 0.) return 0.;   
   else {
     double pi = 3.1415926535;     
     return (-1.*pi)/(exp(lnGamma)*sin(pi*x)*x);
   }   
}
//_____________________________________________________________________________

double Clebsch(const Spin &__j1,const Spin &__m1,const Spin &__j2,
	       const Spin &__m2,const Spin &__J,const Spin &__M){

  if((__m1 + __m2) != __M) return 0.;
  if(abs(__m1) > __j1) return 0;
  if(abs(__m2) > __j2) return 0;

  // convert to pure integers (each 2*spin)
  int j1 = (int)(2.*__j1);
  int m1 = (int)(2.*__m1);
  int j2 = (int)(2.*__j2);
  int m2 = (int)(2.*__m2);
  int J = (int)(2.*__J);
  int M = (int)(2.*__M);

  double n0,n1,n2,n3,n4,n5,d0,d1,d2,d3,d4,A,exp;
  int nu = 0;
  
  double sum = 0;
  while(((d3=(j1-j2-M)/2+nu) < 0)||((n2=(j1-m1)/2+nu) < 0 )) { nu++;}
  while (((d1=(J-j1+j2)/2-nu) >= 0) && ((d2=(J+M)/2-nu) >= 0) 
	 &&((n1=(j2+J+m1)/2-nu) >= 0 )){
    d3=((j1-j2-M)/2+nu);
    n2=((j1-m1)/2+nu);
    d0=dfact((double) nu);
    exp=nu+(j2+m2)/2;
    n0 = (double) pow(-1.,exp);
    sum += ((n0*dfact(n1)*dfact(n2))/(d0*dfact(d1)*dfact(d2)*dfact(d3)));
    nu++;
  }

  if (sum == 0) return 0;

  n0 = J+1;
  n1 = dfact((double) (J+j1-j2)/2);
  n2 = dfact((double) (J-j1+j2)/2);
  n3 = dfact((double) (j1+j2-J)/2);
  n4 = dfact((double) (J+M)/2);
  n5 = dfact((J-M)/2);
  
  d0 = dfact((double) (j1+j2+J)/2+1);
  d1 = dfact((double) (j1-m1)/2);
  d2 = dfact((double) (j1+m1)/2);
  d3 = dfact((double) (j2-m2)/2);
  d4 = dfact((double) (j2+m2)/2);
  
  A = ((double) (n0*n1*n2*n3*n4*n5))/((double) (d0*d1*d2*d3*d4));


  return sqrt(A)*sum;		
}
//_____________________________________________________________________________

double Wigner_d(const Spin &__j,const Spin &__m,const Spin &__n,double __beta){

  int J = (int)(2.*__j);
  int M = (int)(2.*__m);
  int N = (int)(2.*__n);
  int temp_M, k, k_low, k_hi;
  double const_term = 0.0, sum_term = 0.0, d = 1.0;
  int m_p_n, j_p_m, j_p_n, j_m_m, j_m_n;
  int kmn1, kmn2, jmnk, jmk, jnk;
  double kk;

  if (J < 0 || abs (M) > J || abs (N) > J) {
    cerr << endl;
    cerr << "d: you have entered an illegal number for J, M, N." << endl;
    cerr << "Must follow these rules: J >= 0, abs(M) <= J, and abs(N) <= J." 
	 << endl;
    cerr << "J = " << J <<  " M = " << M <<  " N = " << N << endl;
    return 0.;
  }
  
  if (__beta < 0) {
    __beta = fabs (__beta);
    temp_M = M;
    M = N;
    N = temp_M;
  }


  m_p_n = (M + N) / 2;
  j_p_m = (J + M) / 2;
  j_m_m = (J - M) / 2;
  j_p_n = (J + N) / 2;
  j_m_n = (J - N) / 2;
  
  kk = (double)factorial(j_p_m)*(double)factorial(j_m_m)
    *(double)factorial(j_p_n) * (double)factorial(j_m_n) ;
  const_term = pow((-1.0),(j_p_m)) * sqrt(kk);	
  
  k_low = MAX(0, m_p_n);
  k_hi = MIN(j_p_m, j_p_n);

  for (k = k_low; k <= k_hi; k++) {
    
    kmn1 = 2 * k - (M + N) / 2;
    jmnk = J + (M + N) / 2 - 2 * k;
    jmk = (J + M) / 2 - k;
    jnk = (J + N) / 2 - k;
    kmn2 = k - (M + N) / 2;
	
    sum_term += pow ((-1.0), (k)) *
      ((pow (cos (__beta / 2.0), kmn1)) * (pow (sin (__beta / 2.0), jmnk))) /
      (factorial (k) * factorial (jmk) * factorial (jnk) * factorial (kmn2));

  }

  d = const_term * sum_term;
  return d;
}


// double Wigner_d_new(const Spin &__j,const Spin &__m,const Spin &__n,double __beta){

//   int J = (int)(2.*__j);
//   int M = (int)(2.*__m);
//   int N = (int)(2.*__n);
//   int temp_M, k, k_low, k_hi;
//   double const_term = 0.0, sum_term = 0.0, d = 1.0;
//   int m_p_n, j_p_m, j_p_n, j_m_m, j_m_n, m_n_n;
//   int kmn2, jmnk, jmk, jnk;
//   double kk;

//   if (J < 0 || abs (M) > J || abs (N) > J) {
//     cerr << endl;
//     cerr << "d: you have entered an illegal number for J, M, N." << endl;
//     cerr << "Must follow these rules: J >= 0, abs(M) <= J, and abs(N) <= J." 
// 	 << endl;
//    cerr << "J = " << J <<  " M = " << M <<  " N = " << N << endl;
//     return 0.;
//   }
  
//   double preFactor=1.;

//   if (M<0 && N<0){ 
//     preFactor=pow(-1.,M-N);
//     M=-M;
//     N=-N;
//   }

//   if (M < N){
//     __beta = -__beta;
//     temp_M = M;
//     M = N;
//     N = temp_M;
//   }


  
//   m_p_n = (M + N) / 2;
//   j_p_m = (J + M) / 2;
//   j_m_m = (J - M) / 2;
//   j_p_n = (J + N) / 2;
//   j_m_n = (J - N) / 2;
//   m_n_n = (M - N) / 2; 
 
//   kk = (double)factorial(j_p_m)*(double)factorial(j_m_m)
//     *(double)factorial(j_p_n) * (double)factorial(j_m_n) ;
// //   const_term = pow((-1.0),(m_n_n)) * sqrt(kk);
  
//   const_term = sqrt(kk);
  

// //   k_low = MAX(0, 0);
//   k_low = 0;

//   k_hi = MIN(j_p_n, j_m_m);



//   for (k = k_low; k <= k_hi; k++) {
    
//     //    kmn1 = 2 * k - (M - N) / 2;
//     jmnk = J + (N - M) / 2 - 2 * k;
//     jmk = (J - M) / 2 - k;
//     jnk = (J + N) / 2 - k;
//     kmn2 = k + (M - N) / 2;

//     int mmn2k=(M - N)/2  + 2 * k;

// //     std::cout << "kmn1=\t" << kmn1 <<"\n"
// // 	      << "jmnk=\t" << jmnk <<"\n"
// // 	      << "jmk=\t" << jmk <<"\n"
// // 	      << "jnk=\t" << jnk <<"\n"
// // 	      << "kmn2=\t" << kmn2 <<"\n"
// // 	      << "mmn2k=\t" << mmn2k <<"\n"
// // 	      << std::endl;
    	
//     sum_term += pow ((-1.0), (k)) *
//       ((pow (cos (__beta / 2.0), jmnk)) * (pow (sin (__beta / 2.0), mmn2k))) /
//       (factorial (k) * factorial (jmk) * factorial (jnk) * factorial (kmn2));

//   }

//   d = preFactor*const_term * sum_term;
//   return d;
// }
//_____________________________________________________________________________

complex<double> BreitWigner(const Vector4<double> &__p4,double __mass,
				   double __width){
  complex<double> i(0.,1.);
  return __mass*__width/(__mass*__mass - __p4*__p4- i*__mass*__width);
}
//
complex<double> BreitWignerBlattW(const Vector4<double> &__p4, double __massA, 
                                  double __massB, double __mass0, double __width, int __LL){

  if( (__massA < 0.) || (__massB < 0.) ||  __mass0<0.) {
    
    std::cout << "Invalid mass values " << __massA << "\t" << __massB << "\t" << __mass0 <<std::endl;
    assert(0);
  }

   if (__LL==0) return BreitWigner(__p4,__mass0,__width);


  complex<double> i(0.,1.);

  double massAB=__p4.Mass();
   if ( (__massA+__massB) > __mass0) return BreitWigner(__p4,__mass0,__width);
   if ( (__massA+__massB) > massAB) return BreitWigner(__p4,__mass0,__width);

  // break up momentum q0
  double q0=sqrt( (__mass0*__mass0-(__massA+__massB)*(__massA+__massB))
                  * (__mass0*__mass0-(__massA-__massB)*(__massA-__massB)) )
    /(2*__mass0);


  // break up momentum qab
  double qab=sqrt( (massAB*massAB-(__massA+__massB)*(__massA+__massB))
                   * (massAB*massAB-(__massA-__massB)*(__massA-__massB)) )
    /(2*massAB);


  // phase space factor rhoM0
  double  rhoM0=2*q0/__mass0;

  // phase space factor rhoMab
  double  rhoMab=2*qab/massAB;
 
  // R=1.5 GeV-1 for intermediate resonances
  BlattWeisskopf blattWk(__LL, 1.5, q0);

  double blattWqq0=blattWk(qab);

  complex<double>  result=__mass0*__width*blattWqq0/
    (__mass0*__mass0 - massAB*massAB - i*(__mass0*__width* blattWqq0* blattWqq0* rhoMab/rhoM0));

  return result;
}
//

complex<double> FlatteA980(const Vector4<double> &__p4, double __mass0, double g_KK, double g_EtaPi){

  complex<double> i(0.,1.);

  const double mPipm=0.13957018;
  const double mEta=0.547853;
  const double mKpm=0.493677;
  const double mK0=0.497614;
  
  double mAB=__p4.Mass();
  
  if( (mAB < 1e-8)){
    std::cout << "Invalid KK mass  " << mAB << " < 1e-8" <<std::endl;
    assert(0);
  }

  //calculate phase-space factors

  double mPiMmEtaSqr=((mPipm-mEta)/mAB)*((mPipm-mEta)/mAB);
  if (mPiMmEtaSqr>1.) mPiMmEtaSqr=1.;

  double mPiMpEtaSqr=((mPipm+mEta)/mAB)*((mPipm+mEta)/mAB);
  if (mPiMpEtaSqr>1.) mPiMpEtaSqr=1.;

  double mKpmmK0Sqr=((mKpm-mK0)/mAB)*((mKpm-mK0)/mAB);
  if (mKpmmK0Sqr>1.) mKpmmK0Sqr=1.;

  double mKpmpK0Sqr=((mKpm+mK0)/mAB)*((mKpm+mK0)/mAB);
  if (mKpmpK0Sqr>1.) mKpmpK0Sqr=1.;
  
  double rho1=sqrt( (1-mPiMmEtaSqr)* (1-mPiMpEtaSqr) );
  double rho2=sqrt( (1-mKpmmK0Sqr)* (1-mKpmpK0Sqr) );

  complex<double>  result=g_KK*g_KK/( __mass0*__mass0 - mAB*mAB - i * (rho1*g_EtaPi*g_EtaPi+rho2*g_KK*g_KK) );
 
  return result; 
}


complex<double> Flatte(const Vector4<double> &__p4, std::pair<const double, const double>& decPair1, std::pair<const double
		       , const double>& decPair2, double __mass0, double g1, double g2){

  complex<double> i(0.,1.);
  const double m1a=decPair1.first;
  const double m1b=decPair1.second;
  const double m2a=decPair2.first;
  const double m2b=decPair2.second;

  double mAB=__p4.Mass();

  //calculate phase-space factors 

  double ph1Sqr=mAB*mAB/4-m1a*m1b;

  complex<double> gamma11(0.,0.);
  if (ph1Sqr>=0.) gamma11=complex<double>( g1*sqrt(ph1Sqr), 0.);
  else gamma11=complex<double>(0., g1*sqrt(-ph1Sqr));


  double ph2Sqr=mAB*mAB/4-m2a*m2b;

  complex<double> gamma22(0.,0.);
  if (ph2Sqr>=0.) gamma22=complex<double>(g2*sqrt(ph2Sqr), 0.);
  else gamma22=complex<double>(0., g2*sqrt(-ph2Sqr));

  complex<double> gammaLow(0.,0.);
  if( (m1a+m1b) < (m2a+m2b) ) gammaLow=gamma11;
  else gammaLow=gamma22;

  complex<double>  result=__mass0*sqrt(gammaLow*gamma11)/( __mass0*__mass0 - mAB*mAB - i * __mass0 * (gamma11+gamma22) );

  return result;

}

//_____________________________________________________________________________

complex<double> ReggePropagator(double __t,double __s,double __a,double __b,
				const Spin &__spin,int __sig,int __exp_fact){
 
  complex<double> prop,numerator,denominator;
  double alpha = __a*__t + __b; // trajectory
  double pi = 3.1415926535; 
  complex<double> i(0.,1.); 
  numerator = pow(__s,alpha - __spin)*pi*__a;
  numerator *= ((double)__sig + ((double)__exp_fact)*exp(-i*pi*alpha));
  double gamma_arg = alpha + 1. - __spin;
  if(gamma_arg < 0.) denominator = 2*pi/(abs(gamma_arg)*Gamma(abs(gamma_arg)));
  else if(gamma_arg == 0.) denominator = 2*pi;
  else denominator = 2*sin(pi*gamma_arg)*Gamma(gamma_arg);

  prop = numerator/denominator;

  return prop;
}
//_____________________________________________________________________________

Spin GetSpin(const string &__spin){
  
  Spin j;
  string::size_type div_pos = __spin.find("/");
  if(div_pos != string::npos){
    string tmp;
    tmp.assign(__spin,0,div_pos);
    double numer = atof(tmp.c_str());
    tmp.assign(__spin,div_pos + 1,__spin.size() - div_pos - 1);
    double denom = atof(tmp.c_str());
    j = numer/denom;
  }
  else j = atof(__spin.c_str());
  
  return j;
}
//_____________________________________________________________________________

void Wigner_d(const Spin &__jmax,double __beta,
	      map<Spin,map<Spin,map<Spin,double> > > &__d){
  __d.clear();
  Spin jmin,one_half = 1/2.;
  if(__jmax == 0.){
    __d[0][0][0] = 1.;
    return;
  }
  double cb = cos(__beta);
  double sb = sin(__beta);
  // j=1 d's
  map<Spin,map<Spin,double> > d1;
  d1[1][1] = (1+cb)/2.;
  d1[1][0] = -sb/sqrt(2.);
  d1[1][-1] = (1-cb)/2.;
  d1[0][1] = sb/sqrt(2.);
  d1[0][0] = cb;
  d1[0][-1] = -sb/sqrt(2.);
  d1[-1][1] = (1-cb)/2.;
  d1[-1][0] = sb/sqrt(2.);
  d1[-1][-1] = (1+cb)/2.;

  if(__jmax.Denominator() == 1){ // integral spins
    __d[0][0][0] = 1.0;
    if(__jmax == 0.) return;
    __d[1][1][1] = d1[1][1];
    __d[1][1][0] = d1[1][0];
    __d[1][1][-1] = d1[1][-1];
    __d[1][0][1] = d1[0][1];
    __d[1][0][0] = d1[0][0];
    __d[1][0][-1] = d1[0][-1];
    __d[1][-1][1] = d1[-1][1];
    __d[1][-1][0] = d1[-1][0];
    __d[1][-1][-1] = d1[-1][-1];
    if(__jmax == 1.) return;
    jmin = 2.;
  }
  else { // half-integral spins      
    __d[one_half][one_half][one_half] = cos(__beta/2.);
    __d[one_half][one_half][-one_half] = -sin(__beta/2.);      
    __d[one_half][-one_half][one_half] = sin(__beta/2.);
    __d[one_half][-one_half][-one_half] = cos(__beta/2.);
    if(__jmax == one_half) return;
    jmin = 3/2.;
  }

  for(Spin j = jmin; j <= __jmax; j++){
    for(Spin m = -j; m <= j; m++){
      for(Spin n = -j; n <= j; n++){
	double djmn = 0.;
	for(Spin mm = -1; mm <= 1; mm++){
	  for(Spin nn = -1; nn <= 1; nn++){
	    djmn += Clebsch(j-1,m-mm,1,mm,j,m)*Clebsch(j-1,n-nn,1,nn,j,n)
	      *__d[j-1][m-mm][n-nn]*d1[mm][nn];
	  }
	}
	__d[j][m][n] = djmn;
      }
    }
  }  
}
//_____________________________________________________________________________

void Wigner_D(const Spin &__jmax,double __alpha,double __beta,double __gamma,
	      map<Spin,map<Spin,map<Spin,complex<double> > > > &__D){

  complex<double> i(0.,1.);
  map<Spin,map<Spin,map<Spin,double> > > d;
  Wigner_d(__jmax,__beta,d);
  Spin jmin;
  if(d.find(0) != d.end()) jmin = 0;
  else jmin = 1/2.;
  for(Spin j = jmin; j <= __jmax; j++){
    for(Spin m = -j; m <= j; m++){
      for(Spin n = -j; n <= j; n++) 
	__D[j][m][n] = exp(-i*(m*__alpha + n*__gamma))*d[j][m][n];
    }
  }
}
//_____________________________________________________________________________
