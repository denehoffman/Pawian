#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include <boost/multi_array.hpp>
#include "Examples/Tutorial/LineShapes/PiPiSWaveTMatrix.hh"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh" 
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/TMatrixBase.hh"
#include "PwaDynamics/TMatrixRel.hh"
#include "PwaDynamics/TMatrixNonRel.hh"
#include "PwaDynamics/KMatrixPiPiS.hh"
#include "PwaDynamics/KPole.hh"
//#include "PwaDynamics/KPoleNonRel.hh"
#include "PwaDynamics/KMatrixSlowAdlerCorRel.hh"
#include "PwaDynamics/KMatrixNonRel.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"


#include "ErrLogger/ErrLogger.hh"

PiPiSWaveTMatrix::PiPiSWaveTMatrix() :
  _theTFile(0)
  ,_invPiPiMassH1(0)
  ,_invPiPiMassRelH1(0)
  ,_argandH2()
  ,_argandRelH2()
{
  std::string rootFileName="./PiPiSWaveTMatrix.root";
  _theTFile=new TFile(rootFileName.c_str(),"recreate");

  int size=2000;
  double massMin=2.*0.135;
  double massMax=1.9;

  double stepSize=(massMax-massMin)/size;

  _invPiPiMassH1= new TH1F("_invPiPiMassH1","inv #pi^{-} #pi^{+} mass",size+1, massMin, massMax);
  _invPiPiMassH1->SetYTitle("|T_{00}|^{2}");
  _invPiPiMassRelH1= new TH1F("_invPiPiMassRelH1","inv #pi^{-} #pi^{+} mass (rel)",size+1, massMin, massMax);
  _invPiPiMassRelH1->SetYTitle("|#rho_{00}#hat{T}_{00}|^{2}");
  _absT00RelH1= new TH1F("_absT00RelH1","|#rho_{00}#hat{T}_{00}| (rel)",size+1, massMin, massMax);
  _absT00RelH1->SetYTitle("|#rho_{00}#hat{T}_{00}|");
  _absS00RelH1= new TH1F("_absS00RelH1","abs(S00) (rel)",size+1, massMin, massMax);
  _absS00RelH1->SetYTitle("|1-2i#rho_{00}#hat{T}_{00}|");
  _sqrT00RelSigmaPoleH1= new TH1F("_sqrT00RelSigmaPoleH1","|#rho_{00}#hat{T}_{00}|^{2} (rel) #{sigma} pole",size+1, massMin, massMax);
  _pipiPhaseSpaceFactorH2= new TH2F("_pipiPhaseSpaceFactorH2","#pi^{-} #pi^{+} phase space factor",size+1, 0., 2., 1000, 0., 1.);
  _pipipipiPhaseSpaceFactorH2= new TH2F("_pipipipiPhaseSpaceFactorH2","4#pi^{-} phase space factor",size+1, 0., 2., 1000, 0., 1.);
  _argandH2=new TH2F("_argandH2","Argand plot K matrix",301, -1., 1., 301, 0., 1.3);
  _argandH2->SetYTitle("Im(T_{00})");
  _argandH2->SetXTitle("Re(T_{00})");
  _argandRelH2=new TH2F("_argandRelH2","Argand plot K matrix (rel)",301, -1., 1., 301, -1.3, 1.3);
  _argandRelH2->SetYTitle("Im(#rho_{00}#hat{T}_{00})");
  _argandRelH2->SetXTitle("Re(#rho_{00}#hat{T}_{00})");
  _phaseShiftH2=new TH2F("_phaseShiftH2", "phase shift",301, massMin, massMax, 301, 0., 3.1415);
  _phaseShiftRelH2=new TH2F("_phaseShiftRelH2", "phase shift (rel)",301, massMin, massMax, 301, 0., 3.1415);



  boost::shared_ptr<KMatrixSlowAdlerCorRel> theKMatrix(new KMatrixPiPiS());
  boost::shared_ptr<TMatrixBase> theTMatrix(new TMatrixRel(theKMatrix));


  // non relativistic TMatrix
  boost::shared_ptr<KMatrixBase> retrievedKMatrix=theTMatrix->kMatrix();
  vector<boost::shared_ptr<AbsPhaseSpace> > thePhpVecs=retrievedKMatrix->phaseSpaceVec();
  vector<boost::shared_ptr<KPole> > thePoles=retrievedKMatrix->kpoles();

  boost::shared_ptr<KMatrixNonRel> theKMatrixNonRel(new KMatrixNonRel(thePoles, thePhpVecs));
  boost::shared_ptr<TMatrixBase> theTMatrixNonRel(new TMatrixNonRel(theKMatrixNonRel));


  vector<boost::shared_ptr<KPole> >  sigmaPoleVec;
  sigmaPoleVec.push_back(thePoles[0]);
  boost::shared_ptr<array_type_2d> fScatProd=theKMatrix->fScatProd();
  double soScat=theKMatrix->s0Scat(); 
  boost::shared_ptr<KMatrixSlowAdlerCorRel> theKMatrixSigmaPole(new KMatrixSlowAdlerCorRel(sigmaPoleVec, thePhpVecs, fScatProd, soScat));
  boost::shared_ptr<TMatrixBase> theTMatrixSigmaPole(new TMatrixRel(theKMatrixSigmaPole));

  for (double mass=massMin; mass<massMax; mass+=stepSize){
    Vector4<double> mass4Vec(mass, 0.,0.,0.);
    
    theTMatrixNonRel->evalMatrix(mass);
        
    complex<double> currentValLow=(*theTMatrixNonRel)(0,0);
    _invPiPiMassH1->Fill(mass4Vec.M(), norm(currentValLow) );    
    _argandH2->Fill(currentValLow.real(),currentValLow.imag());
    _phaseShiftH2->Fill(mass, atan2(currentValLow.imag(), currentValLow.real())); 


    theTMatrix->evalMatrix(mass);    
    complex<double> currentValLowRel=(*theTMatrix)(0,0);
    complex<double> S00_rel=complex<double>(1.,0.)+2.*complex<double>(0.,1.)*thePhpVecs[0]->factor(mass4Vec.M()).real()*(*theTMatrix)(0,0);

    _invPiPiMassRelH1->Fill(mass4Vec.M(), norm(sqrt( (thePhpVecs[0]->factor(mass4Vec.M())).real() )*currentValLowRel*sqrt( (thePhpVecs[0]->factor(mass4Vec.M())).real() ) ) );
    _absT00RelH1->Fill(mass4Vec.M(), sqrt(norm(sqrt( (thePhpVecs[0]->factor(mass4Vec.M())).real() )*currentValLowRel*sqrt( (thePhpVecs[0]->factor(mass4Vec.M())).real() ) )) );
    _absS00RelH1->Fill(mass4Vec.M(), sqrt(norm(S00_rel))); 
       
    _argandRelH2->Fill( sqrt( (thePhpVecs[0]->factor(mass4Vec.M())).real() )*currentValLowRel.real()*sqrt( (thePhpVecs[0]->factor(mass4Vec.M())).real() ),sqrt( (thePhpVecs[0]->factor(mass4Vec.M())).real() )*currentValLowRel.imag()*sqrt( (thePhpVecs[0]->factor(mass4Vec.M())).real() )); 
    _phaseShiftRelH2->Fill(mass, atan2(currentValLowRel.imag(), currentValLowRel.real()));

    _pipiPhaseSpaceFactorH2->Fill(mass, thePhpVecs[0]->factor(mass4Vec.M()).real() );
    _pipipipiPhaseSpaceFactorH2->Fill(mass, thePhpVecs[2]->factor(mass4Vec.M()).real() );

    theTMatrixSigmaPole->evalMatrix(mass);
    complex<double> currentValRelSigmaPole=(*theTMatrixSigmaPole)(0,0);
    _sqrT00RelSigmaPoleH1->Fill( mass4Vec.M(), norm( sqrt((thePhpVecs[0]->factor(mass4Vec.M())).real())*currentValRelSigmaPole*sqrt( (thePhpVecs[0]->factor(mass4Vec.M())).real())) );     
  }
}

PiPiSWaveTMatrix::~PiPiSWaveTMatrix()
{
  _theTFile->Write();
  _theTFile->Close();
}

