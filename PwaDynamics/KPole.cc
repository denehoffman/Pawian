#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"

KPole::KPole(vector<double>& g_i, double mass_0):
  Matrix< complex<double> >::Matrix(int(g_i.size()), int(g_i.size()))
  , _g_i(g_i)
  , _poleMass(mass_0)
{
}

KPole::~KPole(){
}

void KPole::evalMatrix(const double mass){
  for (int i=0; i< int(_g_i.size()); ++i){
    for (int j=0; j< int(_g_i.size()); ++j){
      this->operator()(i,j)= ( _g_i[i]*_g_i[j])/(_poleMass*_poleMass-mass*mass);
     }
   }
}

