#include "qft++/relativistic-quantum-mechanics/BlattWeisskopf.hh"

#include <iostream>
#include <assert.h>
#include <math.h>
using std::endl;

BlattWeisskopf::BlattWeisskopf(int LL, double R, double p0)
  : _LL(LL), _radial(R), _p0(p0)
{
  if(R < 0) {

    std::cout << "Radius " << R << " negative" << std::endl;
    assert(0);
  }

  _radial = R;

  // compute formula for nominal momentum

  _F0 = compute(_p0);
  if(_F0 <= 0) {
    
    std::cout << "Invalid nominal form factor computed " << _F0 << std::endl;
    assert(0);
  } 
}

BlattWeisskopf::BlattWeisskopf(const BlattWeisskopf& other)
  : _LL(other._LL), _radial(other._radial), _p0(other._p0), _F0(other._F0)
{}

BlattWeisskopf::~BlattWeisskopf()
{}

double BlattWeisskopf::operator()(double p) const
{
  double ret = compute(p)/_F0;
  //  report(INFO,"EvtGen") << p << " " << _p0 << " " << _F0 << " " << _LL << " " << _radial << " " << ret << endl;
  return ret;
}

// Blatt-Weisskopf form factors
// see e.g. hep-ex/0011065
// Dalitz Analysis of the Decay D0->K-pi+pi0 (CLEO)
//
// p   - momentum of either daugher in the meson rest frame,
//       the mass of the meson is used
// pAB - momentum of either daughter in the candidate rest frame
//       the mass of the candidate is used
// R - meson radial parameter
// 
// In the CLEO paper R=5 GeV-1 for D0, R=1.5 for intermediate resonances

double BlattWeisskopf::compute(double p) const
{

  double result=1.;
  if(p < 0) {
    
    std::cout << "Momentum " << p << " negative in form factor calculation" << std::endl;
    assert(0);
  }
  else {
    
    double x = p*p*_radial*_radial;
    
//     if(0 == _LL) return 1.;
//     else
//       if(1 == _LL) return sqrt(1.0/(1.0+x));
//       else
// 	if(2 == _LL) return sqrt(1.0/(1.0+x/3.0+x*x/9.0));
// 	else {
// 	  std::cout << "Angular momentum " << _LL << " not implemented" << std::endl;
// 	  assert(0);
// 	}
//   }

    if(0 == _LL) result=1.;
    else
      if(1 == _LL) result=sqrt(x/(1.0+x));
      else
	if(2 == _LL) result=sqrt((13.*x*x)/((x-3.)*(x-3.)+9.*x));
	else
	  if(3 == _LL) result=sqrt((277.*x*x*x)/(x*(x-15.)*(x-15.)+9.*(2.*x-5)*(2.*x-5)));
	  else {
	    std::cout << "Angular momentum " << _LL << " not implemented" << std::endl;
	    assert(0);
	  }
  }

  return result;
}

