#include "PspGen/EvtGenKine.hh"
#include "EvtVector4R.hh"
#include "EvtRandom.hh"
#include "EvtRandomEngine.hh"
#include "EvtSimpleRandomEngine.hh"

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/WeightContainer.h"
#include "HepMC/Units.h"
#include "HepMC/IO_GenEvent.h"

#include <TROOT.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include <TApplication.h>

#include <iostream>

int main(int argc, char* argv[])
{
  HepMC::IO_GenEvent hepMC_out("HepMCEvt.out",std::ios::out);

  EvtSimpleRandomEngine myRandom(4711);
  EvtRandom::setRandomEngine(&myRandom);

  // variables for the production e+ e- Psi2s
  double mass_e=0.511e-3;
  double mass_psi2s = 3.68609;
  double pe_tot=0.5*sqrt(mass_psi2s*mass_psi2s-2.*mass_e);

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

  for (int count = 0; count < 1000000; count++) {

    EvtGenKine::PhaseSpace(firstNdaug, firstMass, firstP4, firstMp);
    EvtGenKine::PhaseSpace(ndaug, mass, p4, firstP4[0].mass());
    for (int t = 0; t < ndaug; t++) {
      p4[t].applyBoostTo(firstP4[0]);
    }

    //    dalitz.Fill((p4[0]+p4[1])*(p4[0]+p4[1]), (p4[1]+p4[2])*(p4[1]+p4[2]));

     HepMC::GenEvent* evt = new HepMC::GenEvent( 20, count );
    evt->use_units(HepMC::Units::GEV, HepMC::Units::MM);

    // create production vertex
    HepMC::GenVertex* vtx_prod = new HepMC::GenVertex();
    evt->add_vertex( vtx_prod );
    HepMC::GenParticle* eplus_part = new HepMC::GenParticle( HepMC::FourVector(0., 0., pe_tot, sqrt(mass_e*mass_e+pe_tot*pe_tot)), -11, 4 );
    HepMC::GenParticle* eminus_part = new HepMC::GenParticle( HepMC::FourVector(0., 0., -pe_tot, sqrt(mass_e*mass_e+pe_tot*pe_tot)), 11, 4 );
    vtx_prod->add_particle_in( eplus_part );
    vtx_prod->add_particle_in( eminus_part ); 

    HepMC::GenParticle* psi2s_part=new HepMC::GenParticle( HepMC::FourVector(0., 0., 0., mass_psi2s), 20443, 2 );
    vtx_prod->add_particle_out(psi2s_part);

    // create Psi(2S) vertex
    HepMC::GenVertex* vtx_psi2S = new HepMC::GenVertex();
    evt->add_vertex( vtx_psi2S );
    vtx_psi2S->add_particle_in( psi2s_part );
    HepMC::GenParticle* chi_c0_part = new HepMC::GenParticle(HepMC::FourVector(firstP4[0].get(1), firstP4[0].get(2), firstP4[0].get(3), firstP4[0].get(0)), 10441, 2 );
    vtx_psi2S->add_particle_out(chi_c0_part);
    vtx_psi2S->add_particle_out( new HepMC::GenParticle( HepMC::FourVector(firstP4[1].get(1), firstP4[1].get(2), firstP4[1].get(3), firstP4[1].get(0)), 22, 1 ) );


    // now do the chi_c0 vertex
    HepMC::GenVertex* vtx_chi_c0 = new HepMC::GenVertex();
    evt->add_vertex( vtx_chi_c0 );
    vtx_chi_c0->add_particle_in( chi_c0_part );
    vtx_chi_c0->add_particle_out( new HepMC::GenParticle( HepMC::FourVector(p4[0].get(1),p4[0].get(2),p4[0].get(3), p4[0].get(0)), 321, 1 ));
    vtx_chi_c0->add_particle_out( new HepMC::GenParticle( HepMC::FourVector(p4[1].get(1),p4[1].get(2),p4[1].get(3), p4[1].get(0)), -321, 1 ) );
    vtx_chi_c0->add_particle_out( new HepMC::GenParticle( HepMC::FourVector(p4[2].get(1),p4[2].get(2),p4[2].get(3), p4[2].get(0)), 111, 2 ) );
    vtx_chi_c0->add_particle_out( new HepMC::GenParticle( HepMC::FourVector(p4[3].get(1),p4[3].get(2),p4[3].get(3), p4[3].get(0)), 111, 2 ) );

    
    evt->print(std::cout);
    hepMC_out << evt;
    delete evt;

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

//   massKpipi.Draw();
//   rootapp->Run();

  return 0;
}
