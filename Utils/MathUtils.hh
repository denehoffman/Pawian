#ifndef MATHUTILS_HH
#define MATHUTILS_HH

#include <cassert>

inline double pow(double x, int p)
{
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

inline double pow(double x, unsigned int p)
{
  if(p == 0) return 1.0;
  if(x == 0.0) return 0.0;
	
  double r = 1.0;
  for(;;) {
    if(p & 1) r *= x;
    if((p >>= 1) == 0)	return r;
    x *= x;
  }
}


inline int pow(int x, int p)
{
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

inline unsigned int pow(unsigned int x, unsigned int p)
{
  if(p == 0) return 1;
  if(x == 0) return 0;
	
  unsigned int r = 1;
  for(;;) {
    if(p & 1) r *= x;
    if((p >>= 1) == 0)	return r;
    x *= x;
  }
}

#endif
