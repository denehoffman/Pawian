#include "PwaDynamics/PhaseSpace4Pi.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"

PhaseSpace4Pi::PhaseSpace4Pi():
  AbsPhaseSpace()
  , _piMass(0.13957)
{
}

PhaseSpace4Pi::~PhaseSpace4Pi(){

}

complex<double> PhaseSpace4Pi::factor(const double mass){
  double mass_sqr=mass*mass;
  complex<double> result(0.,0.);
  if( mass_sqr <= 1 )
    {
      double real   = 1.2274 + .00370909 / ( mass_sqr * mass_sqr ) - .111203 / mass_sqr - 6.39017 * mass_sqr + 16.8358*mass_sqr*mass_sqr - 21.8845*mass_sqr*mass_sqr*mass_sqr + 11.3153*mass_sqr*mass_sqr*mass_sqr*mass_sqr;
      double cont32 = sqrt(1.0-(16.0*_piMass*_piMass));
      result = complex<double>( cont32 * real, 0 );
    }
  else result = complex<double>( sqrt( 1. - 16. * _piMass * _piMass / mass_sqr ), 0. );
  return result;
}

complex<double> PhaseSpace4Pi::breakUpMom(const double mass){
  return factor(mass)*mass/2.;
}

