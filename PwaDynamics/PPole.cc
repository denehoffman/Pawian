#include "PwaDynamics/PPole.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"

PPole::PPole(complex<double>& beta, vector<double>& g_i, double mass_0):
  KPole(g_i, mass_0, g_i.size(), 1)
  ,_beta(beta)
{
}

PPole::~PPole(){
}

void PPole::evalMatrix(const double mass){

  for (int i=0; i< int(_g_i.size()); ++i){
    this->operator()(i,0)= (_beta*_g_i[i])/(_poleMass*_poleMass-mass*mass);
  }
}

