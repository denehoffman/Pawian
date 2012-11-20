#include "PwaDynamics/KPoleBarrier.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"

KPoleBarrier::KPoleBarrier(vector<double>& g_i, double mass_0, vector<boost::shared_ptr<AbsPhaseSpace> > phpVecs, int orbMom):
  KPole(g_i, mass_0)
  , _phpVecs(phpVecs)
  , _orbMom(orbMom)
{
}

KPoleBarrier::~KPoleBarrier(){
}

void KPoleBarrier::evalMatrix(const double mass){

  vector< complex<double> > barrierFactor;

  for (int i=0; i< int(_phpVecs.size()); ++i){
    complex<double> breakUpM=_phpVecs[i]->breakUpMom(mass);
    complex<double> breakUpM0=_phpVecs[i]->breakUpMom(_poleMass);
    barrierFactor.push_back(breakUpM/breakUpM0);
  }

  for (int i=0; i< int(_g_i.size()); ++i){
    for (int j=0; j< int(_g_i.size()); ++j){
      this->operator()(i,j)= ( _g_i[i]*barrierFactor[i]*_g_i[j]*barrierFactor[j])/(_poleMass*_poleMass-mass*mass);
     }
   }
}

