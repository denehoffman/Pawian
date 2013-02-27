#include "PwaDynamics/Flatte.hh"

#include "Particle/Particle.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"


Flatte::Flatte(std::pair<Particle*, Particle*>& decPair1, std::pair<Particle*, Particle*>& decPair2) :
  _mass11(decPair1.first->mass())
  ,_mass12(decPair1.second->mass())
  ,_mass21(decPair2.first->mass())
  ,_mass22(decPair2.second->mass())
 {
 }

Flatte::Flatte(std::pair <const double, const double>& massPair1, std::pair <const double, const double>& massPair2) :
  _mass11(massPair1.first)
  ,_mass12(massPair1.second)
  ,_mass21(massPair2.first)
  ,_mass22(massPair2.second)
 {
 }
Flatte::~Flatte(){
}

complex<double> Flatte::calcFirstChannel(double currentMass, double mass0, double g1, double g2){
  
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

complex<double> Flatte::calcSecondChannel(double currentMass, double mass0, double g1, double g2){
  
  complex<double> i(0.,1.);
  
  //calculate gammas with phase-space factors 
  complex<double> gamma11=g1*breakupMomQ(currentMass, _mass11, _mass12);
  complex<double> gamma22=g2*breakupMomQ(currentMass, _mass21, _mass22);

  complex<double> gammaLow(0.,0.);
  if( (_mass11+_mass12) < (_mass21+_mass22) ) gammaLow=gamma11;
  else gammaLow=gamma22;

  complex<double>  result=mass0*sqrt(gammaLow*gamma22)/( mass0*mass0 - currentMass*currentMass - i * mass0 * (gamma11+gamma22) );
  return result;
}
