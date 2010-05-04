#include "Particle/Parity.hh"

Parity::Parity() :
  pvalue(-2)
{
}

Parity::Parity(int value) 
{
  if (value > -2 && value < 2)
    pvalue = value;
  else
    pvalue = -2;
}

Parity::~Parity()
{
}

int Parity::parity()
{
  return pvalue;
}
