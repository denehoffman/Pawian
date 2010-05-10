#ifndef ERRVALUE_HH
#define ERRVALUE_HH

#include <iostream>

class ErrValue
{
public:

  ErrValue();
  ErrValue(double value, double errPlus = 0.0, double errMinus = 0.0);
  ~ErrValue();
  
  double mean() const;
  double plusErr() const;
  double minusErr() const;

  void mean(double v);
  void plusErr(double v);
  void minusErr(double v);

  ErrValue operator+(ErrValue& v) const;
  ErrValue operator-(ErrValue& v) const;
  ErrValue operator*(ErrValue& v) const;
  ErrValue operator/(ErrValue& v) const;

  ErrValue& operator=(ErrValue& v);
  ErrValue& operator+=(ErrValue& v);
  ErrValue& operator-=(ErrValue& v);
  ErrValue& operator*=(ErrValue& v);
  ErrValue& operator/=(ErrValue& v);

  bool operator<(ErrValue& v) const;
  bool operator<=(ErrValue& v) const;
  bool operator>(ErrValue& v) const;
  bool operator>=(ErrValue& v) const;

private:

  double meanValue;
  double plus;
  double minus;

friend std::ostream &operator<<(std::ostream &o, const ErrValue &v);

};

#endif
