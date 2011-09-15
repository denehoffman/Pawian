#include "PspGen/EvtGenKine.hh"
#include "EvtVector4R.hh"
#include "EvtRandom.hh"
#include "EvtRandomEngine.hh"
#include "EvtSimpleRandomEngine.hh"

#include <TROOT.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include <TApplication.h>

#include <iostream>

int main(int argc, char* argv[])
{
  EvtSimpleRandomEngine myRandom;
  EvtRandom::setRandomEngine((EvtRandomEngine*)(&myRandom));

  // variables for the first decay Psi2s->Chi_c0+gamma
  int firstNdaug = 2;
  double firstMass[8];
  EvtVector4R firstP4[8]; 
  double firstMp = 3.68609;

  firstMass[0] = 3.41475;
  firstMass[10] = 0.;

  // variables for the second decay Chi->c0->K+ K- pi0 pi0
  int ndaug = 4;
  double mass[30];
  EvtVector4R p4[30]; 

  mass[0] = .493677;
  mass[1] = .493677;
  mass[2] = .1349766;
  mass[3] = .1349766;

  // 5 events for visual control
  for (int count = 0; count < 5; count++) {
    std::cout << "Event #" << count << std::endl;
    EvtGenKine::PhaseSpace(firstNdaug, firstMass, firstP4, firstMp);
    for (int t = 0; t < firstNdaug; t++) {
      std::cout << firstP4[t] << " m = " << sqrt(firstP4[t]*firstP4[t]) << std::endl;
    }
    // now use the mass of the Chi_c0 from the first decay as input for the second decay
    EvtGenKine::PhaseSpace(ndaug, mass, p4, firstP4[0].mass());
    for (int t = 0; t < ndaug; t++) {
      std::cout << p4[t] << " m = " << sqrt(p4[t]*p4[t]) << std::endl;
      p4[t].applyBoostTo(firstP4[0]);
      std::cout << p4[t] << " boostet: m = " << sqrt(p4[t]*p4[t]) << std::endl;
    }
  }

  TApplication* rootapp = new TApplication("example",&argc, argv);
  TH2F dalitz("Dalitz plot", "Dalitz Plot", 80,0.,5., 80,0.,6.);

  TH1F massKpi("invMassKpi", "invMassKpi", 512, 0., 3.);
  TH1F massKpipi("invMassKpipi", "invMassKpipi", 512, 0., 3.);
  TH1F cosThetapi("cosThetapi", "cosThetapi", 512, -1.2, 1.2);
  TH1F cosThetaKpi("cosThetaKpi", "cosThetaKpi", 512, -1.2, 1.2);

  for (int count = 0; count < 5000000; count++) {

    EvtGenKine::PhaseSpace(firstNdaug, firstMass, firstP4, firstMp);
    EvtGenKine::PhaseSpace(ndaug, mass, p4, firstP4[0].mass());
    for (int t = 0; t < ndaug; t++) {
      p4[t].applyBoostTo(firstP4[0]);
    }

    //    dalitz.Fill((p4[0]+p4[1])*(p4[0]+p4[1]), (p4[1]+p4[2])*(p4[1]+p4[2]));

    massKpi.Fill( (p4[0]+p4[2]).mass(), 1.);
    massKpi.Fill( (p4[0]+p4[3]).mass(), 1.);
    massKpi.Fill( (p4[1]+p4[2]).mass(), 1.);
    massKpi.Fill( (p4[1]+p4[3]).mass(), 1.);

    massKpipi.Fill( (p4[0]+p4[2]+p4[3]).mass(), 1.);
    massKpipi.Fill( (p4[1]+p4[2]+p4[3]).mass(), 1.);

    cosThetapi.Fill( p4[2].get(3) /
		     sqrt(p4[2].get(1)*p4[2].get(1) + p4[2].get(2)*p4[2].get(2) + p4[2].get(3)*p4[2].get(3)), 
		     1.);

    cosThetaKpi.Fill( (p4[0]+p4[2]).get(3) /
		     sqrt((p4[0]+p4[2]).get(1)*(p4[0]+p4[2]).get(1) + 
			  (p4[0]+p4[2]).get(2)*(p4[0]+p4[2]).get(2) + 
			  (p4[0]+p4[2]).get(3)*(p4[0]+p4[2]).get(3)), 
		     1.);
  }

  massKpipi.Draw();
  rootapp->Run();

  return 0;
}
