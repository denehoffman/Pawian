#include "PwaDynamics/FlatteShape.hh"

#include "Particle/Particle.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"


FlatteShape::FlatteShape(std::pair<Particle*, Particle*>& decPair1, std::pair<Particle*, Particle*>& decPair2) :
  _mass11(decPair1.first->mass())
  ,_mass12(decPair1.second->mass())
  ,_mass21(decPair2.first->mass())
  ,_mass22(decPair2.second->mass())
 {
 }

FlatteShape::~FlatteShape(){
}

complex<double> FlatteShape::calculate(double currentMass, double mass0, double g1, double g2){
  
  complex<double> i(0.,1.);
  
  //calculate gammas with phase-space factors 
  complex<double> gamma11=g1*breakupMomQ(currentMass, _mass11, _mass12);
  complex<double> gamma22=g2*breakupMomQ(currentMass, _mass21, _mass22);

  complex<double> gammaLow(0.,0.);
  if( (_mass11+_mass12) < (_mass21+_mass22) ) gammaLow=gamma11;
  else gammaLow=gamma22;

  complex<double>  result=mass0*sqrt(gammaLow*gamma11)/( mass0*mass0 - currentMass*currentMass - i * mass0 * (gamma11+gamma22) );
  return result;
}

