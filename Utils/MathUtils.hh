#pragma once
#include <iostream>
#include <cassert>
#include <math.h>

namespace MathUtils{
   
    inline double Gauss(double __x, double __mean, double __sigma) {
        __sigma= std::abs(__sigma) ;
        double arg = (__x-__mean)/__sigma;
        double res = exp(-0.5*arg*arg);
        return res/(2.50662827463100024*__sigma); 
    }

  inline int factorial(int __i) {
    if (__i<0){
      std::cerr << std::endl;
      std::cerr << "factoral value " << __i << " must be >=0 !!! "  
		<< std::endl;
      
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

  inline int double_factorial(int n){
    if (n<-1){
      std::cerr << std::endl;
      std::cerr << "factoral value " << n << " must be >=-1 !!! "
                << std::endl;
      
      exit(1);
    }
    
    if (n==-1) return 1;
    int result = 1;
    
    if ( n % 2 == 0){ 
      for(int i =1; i<=floor(n/2.); i++){
	result*=2*i;
      }
    }
    
    if ( n % 2 == 1){ 
      for(int i =1; i<=floor((n+1.)/2.); i++){
	result*=2*i-1;
      }
    }
            
            return result;
  }

  inline double pow(double x, int p) {
    if(p == 0) return 1.0;
    if(x == 0.0 && p > 0) return 0.0;
    if(p < 0) {p=-p; x=1/x;}
    
    double r = 1.0;
    for(;;) {
      if(p & 1) r *= x;
      if((p >>= 1) == 0)	return r;
      x *= x;
    }
  }
  
  inline double pow(double x, unsigned int p) {
    if(p == 0) return 1.0;
    if(x == 0.0) return 0.0;
    
    double r = 1.0;
    for(;;) {
      if(p & 1) r *= x;
      if((p >>= 1) == 0)	return r;
      x *= x;
    }
  }
  
  inline int pow(int x, int p) {
    if(p == 0) return 1;
    if(x == 0 && p > 0) return 0;
    if(p < 0) {assert(x == 1 || x == -1); return (-p % 2) ? x : 1;}
    
    int r = 1;
    for(;;) {
      if(p & 1) r *= x;
      if((p >>= 1) == 0)	return r;
      x *= x;
    }
  }
  
  inline unsigned int pow(unsigned int x, unsigned int p) {
    if(p == 0) return 1;
    if(x == 0) return 0;
    
    unsigned int r = 1;
    for(;;) {
      if(p & 1) r *= x;
      if((p >>= 1) == 0)	return r;
      x *= x;
    }
  }
  
}
