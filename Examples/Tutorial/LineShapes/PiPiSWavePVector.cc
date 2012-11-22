#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include <boost/multi_array.hpp>
#include "Examples/Tutorial/LineShapes/PiPiSWavePVector.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/TMatrixBase.hh"
#include "PwaDynamics/TMatrixRel.hh"
#include "PwaDynamics/TMatrixNonRel.hh"
#include "PwaDynamics/KPole.hh"
//#include "PwaDynamics/KPoleNonRel.hh"
#include "PwaDynamics/KMatrixSlowAdlerCorRel.hh"
#include "PwaDynamics/KMatrixNonRel.hh"

#include "PwaDynamics/PPole.hh"
#include "PwaDynamics/PVectorRel.hh"
#include "PwaDynamics/PVectorSlowCorRel.hh"
#include "PwaDynamics/FVector.hh"
#include "PwaDynamics/FVectorPiPiS.hh"


#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"


#include "ErrLogger/ErrLogger.hh"

PiPiSWavePVector::PiPiSWavePVector() :
  _theTFile(0)
{
  std::string rootFileName="./PiPiSWavePVector.root";
  _theTFile=new TFile(rootFileName.c_str(),"recreate");

  std::vector< complex<double> > fProd;
  fProd.resize(5);
  fProd[0]=complex <double>(- 2.19449   , -7.62666);
  fProd[1]=complex <double>(1.87981, -.628378)*fProd[0];
  fProd[2]=complex <double>(4.3242 , 2.75019)*fProd[0];
  fProd[3]=complex <double>(3.22336,  .271048)*fProd[0];
  fProd[4]=complex <double>(0.,0.);

  std::vector<complex <double> > betaPars;
  betaPars.push_back(complex <double>(1.8223,-9.11972));
  betaPars.push_back(complex <double>(-10.1735,-3.88488));
  betaPars.push_back(complex <double>(-23.6712,5.03146));
  betaPars.push_back(complex <double>(-0.0799351,9.15965));
  betaPars.push_back(complex <double>(0.,0.));

  int size=2000;
  double massMin=2.*0.135;
  double massMax=1.9;

  // double massMax=1.4;

  double stepSize=(massMax-massMin)/size;

  _absF0H1= new TH1F("_absF0H1","|#rho_{00}#hat{F}_{0}|",size+1, massMin, massMax);
  _absF0H1->SetYTitle("|#rho_{00}#hat{F}_{0}|");
  _sqrF0H1= new TH1F("_sqrF0H1","|#rho_{00}#hat{F}_{0}|^{2}",size+1, massMin, massMax);
  _sqrF0H1->SetYTitle("|#rho_{00}#hat{F}_{0}|^{2}");
  _argandH2=new TH2F("_argandH2","Argand plot",200, -10., 10., 200, -10., 10.);
  _argandH2->SetYTitle("Im(F_{0})");
  _argandH2->SetXTitle("Re(F_{0})");
  _phaseShiftH2=new TH2F("_phaseShiftH2", "phase shift",301, massMin, massMax, 301, 0., 3.1415);


  // //  double s0Prod=-3.;
  // boost::shared_ptr<PVectorRel> thePVector(new PVectorSlowCorRel(pPoles, phpVecs, fProdMatr, s0Prod));

  // //boost::shared_ptr<PVectorRel> thePVector(new PVectorRel(pPoles, phpVecs));

  boost::shared_ptr<FVector> theFVector(new FVectorPiPiS());

  std::vector<complex <double> >::const_iterator itBeta;
  for( int i=0; i< int(betaPars.size()); ++i){
    theFVector->updateBeta(i, betaPars[i]);    
  } 

  std::vector<complex <double> >::const_iterator itFprod;
  for( int i=0; i< int(fProd.size()); ++i){
    theFVector->updateFprod(i, fProd[i]);    
  } 

  double s0Prod=-0.0737;
  theFVector->updateS0prod(s0Prod);
  
  for (double mass=massMin; mass<massMax; mass+=stepSize){
    Vector4<double> mass4Vec(mass, 0.,0.,0.);
    
    theFVector->evalMatrix(mass);
    complex<double> currentValFvec=(*theFVector)(0,0);
    _absF0H1->Fill(mass4Vec.M(), sqrt(norm(currentValFvec)));
    _sqrF0H1->Fill(mass4Vec.M(), norm(currentValFvec) );
    _argandH2->Fill( currentValFvec.real(), currentValFvec.imag() );
    _phaseShiftH2->Fill(mass, atan2(currentValFvec.imag(), currentValFvec.real()));
  }
}

PiPiSWavePVector::~PiPiSWavePVector()
{
  _theTFile->Write();
  _theTFile->Close();
}

