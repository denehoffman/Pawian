#include "PwaDynamics/PhaseSpaceIsobar.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"

PhaseSpaceIsobar::PhaseSpaceIsobar(double mass1, double mass2):
  AbsPhaseSpace()
  , _mass1(mass1)
  ,_mass2(mass2)
{
}

PhaseSpaceIsobar::~PhaseSpaceIsobar(){

}

complex<double> PhaseSpaceIsobar::factor(const double mass){
  return phaseSpaceFac(mass,_mass1, _mass2);
}

complex<double> PhaseSpaceIsobar::breakUpMom(const double mass){
  return breakupMomQ(mass,_mass1, _mass2);
}

