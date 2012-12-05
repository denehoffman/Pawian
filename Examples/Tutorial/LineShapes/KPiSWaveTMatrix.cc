#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include <boost/multi_array.hpp>
#include "Examples/Tutorial/LineShapes/KPiSWaveTMatrix.hh"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh" 
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/TMatrixBase.hh"
#include "PwaDynamics/TMatrixRel.hh"
#include "PwaDynamics/TMatrixNonRel.hh"
#include "PwaDynamics/KMatrixBase.hh"
#include "PwaDynamics/KMatrixKPiSFocus.hh"
#include "PwaDynamics/KPole.hh"
//#include "PwaDynamics/KPoleNonRel.hh"
#include "PwaDynamics/KMatrixRel.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"


#include "ErrLogger/ErrLogger.hh"

KPiSWaveTMatrix::KPiSWaveTMatrix() :
  _theTFile(0)
{
  std::string rootFileName="./KPiSWaveTMatrix.root";
  _theTFile=new TFile(rootFileName.c_str(),"recreate");

  int size=2000;
  double massMin=0.135+0.5;
  double massMax=1.8;

  double stepSize=(massMax-massMin)/size;

  _KPiAmpRealH1= new TH1F("_KPiAmpRealH1","K #pi amp Re",size+1, massMin, massMax);
  _KPiAmpRealH1->SetYTitle("K #pi amp Re");
  _KPiAmpImagH1= new TH1F("_KPiAmpImagH1","K #pi amp Im",size+1, massMin, massMax);
  _KPiAmpImagH1->SetYTitle("K #pi amp Im");


  boost::shared_ptr<KMatrixBase> theKMatrixIso12(new KMatrixKPiSFocus(1));
  boost::shared_ptr<KMatrixBase> theKMatrixIso32(new KMatrixKPiSFocus(3));
  // boost::shared_ptr<TMatrixRel> theTMatrix(new TMatrixRel(theKMatrix));
  boost::shared_ptr<TMatrixRel> theTMatrix12(new TMatrixRel(theKMatrixIso12));
  boost::shared_ptr<TMatrixRel> theTMatrix32(new TMatrixRel(theKMatrixIso32));

  vector<boost::shared_ptr<AbsPhaseSpace> > thePhpVecs=theKMatrixIso12->phaseSpaceVec();

  
  for (double mass=massMin; mass<massMax; mass+=stepSize){
    Vector4<double> mass4Vec(mass, 0.,0.,0.);
    
    theTMatrix12->evalMatrix(mass);
    theTMatrix32->evalMatrix(mass);        

    complex<double> currentVal=(*theTMatrix12)(0,0)+0.5*(*theTMatrix32)(0,0);
    _KPiAmpRealH1->Fill(mass, sqrt(thePhpVecs[0]->factor(mass4Vec.M()).real())*currentVal.real());    
    _KPiAmpImagH1->Fill(mass, sqrt(thePhpVecs[0]->factor(mass4Vec.M()).real())*currentVal.imag());   
  }
}

KPiSWaveTMatrix::~KPiSWaveTMatrix()
{
  _theTFile->Write();
  _theTFile->Close();
}

