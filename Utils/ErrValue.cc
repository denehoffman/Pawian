//
// class to hold a double value with asymmetric errors
//

#include "Utils/ErrValue.hh"
#include <math.h>


#define sqr(x) (x*x)

ErrValue::ErrValue():
  meanValue(0.0),
  plus(0.0),
  minus(0.0)
{}

ErrValue::ErrValue(double value, double errPlus, double errMinus):
  meanValue(value),
  plus(errPlus),
  minus(errMinus)
{}

ErrValue::~ErrValue()
{}
  
double ErrValue::mean() const
{
  return meanValue;
}

double ErrValue::plusErr() const
{
  return plus;
}

double ErrValue::minusErr() const
{
  return minus;
}

void ErrValue::mean(double v)
{
  meanValue = v;
  return;
}

void ErrValue::plusErr(double v)
{
  plus = v;
  return;
}
void ErrValue::minusErr(double v)
{
  minus = v;
  return;
}

ErrValue ErrValue::operator+(ErrValue& v) const
{
  ErrValue result(meanValue+v.mean(), 
		  sqrt(plus*plus + v.plusErr()*v.plusErr()), 
		  sqrt(minus*minus + v.minusErr()*v.minusErr()));
  return result;
}

ErrValue ErrValue::operator-(ErrValue& v) const
{
  ErrValue result(meanValue-v.mean(), 
		  sqrt(plus*plus - v.plusErr()*v.plusErr()), 
		  sqrt(minus*minus - v.minusErr()*v.minusErr()));
  return result;
}

ErrValue ErrValue::operator*(ErrValue& v) const
{
  ErrValue result(meanValue*v.mean(), 
		  sqrt(sqr(v.mean()*plus) + sqr(meanValue*v.plusErr())), 
		  sqrt(sqr(v.mean()*minus) + sqr(meanValue*v.minusErr())));
  return result;
}

ErrValue ErrValue::operator/(ErrValue& v) const
{
  ErrValue result(meanValue/v.mean(), 
		  sqrt(sqr(v.mean()*plus) + sqr(meanValue*v.plusErr())), 
		  sqrt(sqr(v.mean()*minus) + sqr(meanValue*v.minusErr())));
  return result;
}

ErrValue& ErrValue::operator=(ErrValue& v)
{
  if ( this == &v )  return *this;

  meanValue = v.mean();
  plus = v.plusErr();
  minus = v.minusErr();
  
  return *this;
}

ErrValue& ErrValue::operator+=(ErrValue& v)
{
  meanValue += v.mean();
  plus = sqrt(plus*plus + v.plusErr()*v.plusErr());
  minus = sqrt(minus*minus + v.minusErr()*v.minusErr());
  
  return *this;
}

ErrValue& ErrValue::operator-=(ErrValue& v)
{
  meanValue -= v.mean();
  plus = sqrt(plus*plus + v.plusErr()*v.plusErr()); 
  minus = sqrt(minus*minus + v.minusErr()*v.minusErr());
  
  return *this;
}

ErrValue& ErrValue::operator*=(ErrValue& v)
{
  meanValue *= v.mean();
  plus = sqrt(sqr(v.mean()*plus) + sqr(meanValue*v.plusErr()));
  minus = sqrt(sqr(v.mean()*minus) + sqr(meanValue*v.minusErr()));

  return *this;
}

ErrValue& ErrValue::operator/=(ErrValue& v)
{
  meanValue /= v.mean();
  plus = sqrt(sqr(v.mean()*plus) + sqr(meanValue*v.plusErr())); 
  minus = sqrt(sqr(v.mean()*minus) + sqr(meanValue*v.minusErr()));

  return *this;
}

bool ErrValue::operator<(ErrValue& v) const
{
  return (meanValue < v.mean());
}

bool ErrValue::operator<=(ErrValue& v) const
{
  return (meanValue <= v.mean());
}

bool ErrValue::operator>(ErrValue& v) const
{
  return (meanValue > v.mean());
}

bool ErrValue::operator>=(ErrValue& v) const
{
  return (meanValue >= v.mean());
}


std::ostream &operator<<(std::ostream &o, const ErrValue &v)
{
  o << v.mean() << "+" << v.plusErr() << "-" << v.minusErr();
  return o;
}
