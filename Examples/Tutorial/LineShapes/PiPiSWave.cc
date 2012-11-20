#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include <boost/multi_array.hpp>
#include "Examples/Tutorial/LineShapes/PiPiSWave.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/TMatrixBase.hh"
#include "PwaDynamics/TMatrixRel.hh"
#include "PwaDynamics/TMatrixNonRel.hh"
#include "PwaDynamics/KPole.hh"
//#include "PwaDynamics/KPoleNonRel.hh"
#include "PwaDynamics/KMatrixSlowAdlerCorRel.hh"
#include "PwaDynamics/KMatrixNonRel.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"


#include "ErrLogger/ErrLogger.hh"

PiPiSWave::PiPiSWave(std::map< int, std::vector<double> >& gFactors, std::vector<double >& poleMasses, std::vector<  boost::shared_ptr<AbsPhaseSpace> > phpVecs) :
  _theTFile(0)
  ,_invPiPiMassH1(0)
  ,_invPiPiMassRelH1(0)
  ,_argandH2()
  ,_argandRelH2()
{
  std::string rootFileName="./PiPiSWave.root";
  _theTFile=new TFile(rootFileName.c_str(),"recreate");

  boost::shared_ptr< boost::multi_array<double, 2> > fScatMatr(new boost::multi_array<double, 2>(boost::extents[gFactors.size()][gFactors.size()]));
  (*fScatMatr)[0][0]=0.23399;
  (*fScatMatr)[0][1]=0.15044;
  (*fScatMatr)[0][2]=-0.20545;
  (*fScatMatr)[0][3]=0.32825;
  (*fScatMatr)[0][4]=0.35412;

  for (int i=1; i<int(gFactors.size()); ++i){
    for (int j=0; j<int(gFactors.size()); ++j){
      (*fScatMatr)[i][j]=0.;
    }
  }
  (*fScatMatr)[1][0]=(*fScatMatr)[0][1];
  (*fScatMatr)[2][0]=(*fScatMatr)[0][2];
  (*fScatMatr)[3][0]=(*fScatMatr)[0][3];
  (*fScatMatr)[4][0]=(*fScatMatr)[0][4];
  double s0Scat=-3.92637;

  int size=2000;
  double massMin=2.*0.135;
  double massMax=1.9;
  // double massMax=1.0;

  double stepSize=(massMax-massMin)/size;

  _invPiPiMassH1= new TH1F("_invPiPiMassH1","inv #pi^{-} #pi^{+} mass",size+1, massMin, massMax);
  _invPiPiMassH1->SetYTitle("|T_{00}|^{2}");
  _invPiPiMassRelH1= new TH1F("_invPiPiMassRelH1","inv #pi^{-} #pi^{+} mass (rel)",size+1, massMin, massMax);
  _invPiPiMassRelH1->SetYTitle("|#rho_{00}#hat{T}_{00}|^{2}");
  _argandH2=new TH2F("_argandH2","Argand plot K matrix",301, -1., 1., 301, 0., 1.3);
  _argandH2->SetYTitle("Im(T_{00})");
  _argandH2->SetXTitle("Re(T_{00})");
  _argandRelH2=new TH2F("_argandRelH2","Argand plot K matrix (rel)",301, -1., 1., 301, -1.3, 1.3);
  _argandRelH2->SetYTitle("Im(#rho_{00}#hat{T}_{00})");
  _argandRelH2->SetXTitle("Re(#rho_{00}#hat{T}_{00})");
  _phaseShiftH2=new TH2F("_phaseShiftH2", "phase shift",301, massMin, massMax, 301, 0., 3.1415);
  _phaseShiftRelH2=new TH2F("_phaseShiftRelH2", "phase shift",301, massMin, massMax, 301, 0., 3.1415);

  std::vector< boost::shared_ptr<KPole> > kPoles;
  // std::vector< boost::shared_ptr<KPole> > kPolesNonRel;
  std::map<int, std::vector<double> >::iterator itgFac;
  for (itgFac=gFactors.begin(); itgFac!=gFactors.end(); ++itgFac){

    boost::shared_ptr<KPole> currentPole(new KPole(itgFac->second, poleMasses[itgFac->first]));
    kPoles.push_back(currentPole);
  }  

  boost::shared_ptr<KMatrixBase> theKMatrixNonRel(new KMatrixNonRel(kPoles, phpVecs));
  boost::shared_ptr<TMatrixBase> theTMatrixNonRel(new TMatrixNonRel(theKMatrixNonRel));

  boost::shared_ptr<KMatrixBase> theKMatrix(new KMatrixSlowAdlerCorRel(kPoles, phpVecs, fScatMatr, s0Scat));
  boost::shared_ptr<TMatrixBase> theTMatrix(new TMatrixRel(theKMatrix));  
  for (double mass=massMin; mass<massMax; mass+=stepSize){
    Vector4<double> mass4Vec(mass, 0.,0.,0.);
    
    theTMatrixNonRel->evalMatrix(mass);
        
    complex<double> currentValLow=(*theTMatrixNonRel)(0,0);
    _invPiPiMassH1->Fill(mass4Vec.M(), norm(currentValLow) );    
    _argandH2->Fill(currentValLow.real(),currentValLow.imag());
    _phaseShiftH2->Fill(mass, atan2(currentValLow.imag(), currentValLow.real())); 


    theTMatrix->evalMatrix(mass);
    
    complex<double> currentValLowRel=(*theTMatrix)(0,0);

    _invPiPiMassRelH1->Fill(mass4Vec.M(), norm(sqrt( (phpVecs[0]->factor(mass4Vec.M())).real() )*currentValLowRel)*sqrt( (phpVecs[0]->factor(mass4Vec.M())).real() ) );        
    _argandRelH2->Fill( sqrt( (phpVecs[0]->factor(mass4Vec.M())).real() )*currentValLowRel.real()*sqrt( (phpVecs[0]->factor(mass4Vec.M())).real() ),sqrt( (phpVecs[0]->factor(mass4Vec.M())).real() )*(phpVecs[0]->factor(mass4Vec.M())*currentValLowRel).imag()*sqrt( (phpVecs[0]->factor(mass4Vec.M())).real() )); 
    _phaseShiftRelH2->Fill(mass, atan2(currentValLowRel.imag(), currentValLowRel.real()));    
  }
   
}

PiPiSWave::~PiPiSWave()
{
  _theTFile->Write();
  _theTFile->Close();
}

