
#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>
#include "PwaDynamics/PhaseSpaceIsobar.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpace4Pi.hh"
#include "Examples/Tutorial/LineShapes/PiPiSWave.hh"
#include "ErrLogger/ErrLogger.hh"

int main(int __argc,char *__argv[]){
  ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);

  const double piMass=0.1349766;
  const double KplusMass=0.493677;
  const double K0Mass= 0.497614;
  const double etaMass= 0.547853;
  const double etaprimeMass= 0.95778;

  std::map<int, std::vector<double> > gFactorsMap;
  std::vector<double> poleMasses;
  poleMasses.resize(5);
  poleMasses[0]=0.651;
  poleMasses[1]=1.2036;
  poleMasses[2]=1.55817;
  poleMasses[3]=1.21;
  poleMasses[4]=1.82206;

  gFactorsMap[0].push_back(0.22889); //pipi
  gFactorsMap[0].push_back(-0.55377); //KK
  gFactorsMap[0].push_back(0.); //4pi
  gFactorsMap[0].push_back(-0.39899); //eta eta
  gFactorsMap[0].push_back(-0.34639); //eta eta'

  gFactorsMap[1].push_back(0.94128); //pipi
  gFactorsMap[1].push_back(0.55095); //KK
  gFactorsMap[1].push_back(0.); //4pi
  gFactorsMap[1].push_back(0.39065); //eta eta
  gFactorsMap[1].push_back(0.31503); //eta eta'

  gFactorsMap[2].push_back(0.36856); //pipi
  gFactorsMap[2].push_back(0.23888); //KK
  gFactorsMap[2].push_back(0.55639); //4pi
  gFactorsMap[2].push_back(0.18340); //eta eta
  gFactorsMap[2].push_back(0.18681); //eta eta'

  gFactorsMap[3].push_back(0.3365); //pipi
  gFactorsMap[3].push_back(0.40907); //KK
  gFactorsMap[3].push_back(0.85679); 
  gFactorsMap[3].push_back(0.19906); //eta eta
  gFactorsMap[3].push_back(-0.00984); //eta eta'

  gFactorsMap[4].push_back(0.18171); //pipi
  gFactorsMap[4].push_back(-0.17558); //KK
  gFactorsMap[4].push_back(-0.79658); //4pi
  gFactorsMap[4].push_back(-0.00355); //eta eta
  gFactorsMap[4].push_back(0.22358); //eta eta'

  boost::shared_ptr<AbsPhaseSpace> pipiPhp(new PhaseSpaceIsobar(piMass, piMass));
  boost::shared_ptr<AbsPhaseSpace> kkPhp(new PhaseSpaceIsobar(KplusMass, K0Mass));
  boost::shared_ptr<AbsPhaseSpace> pipipipiPhp(new PhaseSpace4Pi());
  boost::shared_ptr<AbsPhaseSpace> etaetaPhp(new PhaseSpaceIsobar(etaMass, etaMass));
  boost::shared_ptr<AbsPhaseSpace> etaetapPhp(new PhaseSpaceIsobar(etaMass, etaprimeMass));

  std::vector<  boost::shared_ptr<AbsPhaseSpace> > phpVecs;
  phpVecs.push_back(pipiPhp);
  phpVecs.push_back(kkPhp);
  phpVecs.push_back(pipipipiPhp);
  phpVecs.push_back(etaetaPhp);
  phpVecs.push_back(etaetapPhp);

  PiPiSWave pipiSWave(gFactorsMap, poleMasses, phpVecs);

  return 0;
}

