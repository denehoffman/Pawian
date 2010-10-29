#ifndef BLATT_WEISSKOPF_HH
#define BLATT_WEISSKOPF_HH

class BlattWeisskopf {


public:

  BlattWeisskopf(int LL, double R, double p0);
  BlattWeisskopf(const BlattWeisskopf&);
  ~BlattWeisskopf();

  double operator()(double p) const;

private:

  int    _LL;   // angular momentum of daughters
  double _radial;    // resonance radial parameter
  double _p0;

  double _F0;   // formula evaluated at _p0
  double compute(double p) const;

};

#endif


