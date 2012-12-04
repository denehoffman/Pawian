#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Examples/Tutorial/DfuncClebschG/DfuncHist.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

#include "ErrLogger/ErrLogger.hh"

DfuncHist::DfuncHist(Spin J, Spin lam1, Spin lam2) :
  _J(J),
  _lam1(lam1),
  _lam2(lam2)

{
  if ( fabs(J)< fabs(lam1) || fabs(J) < fabs(lam2) ){
    Alert << " fabs(J)< fabs(lam1) or fabs(J) < fabs(lam2) " << endmsg;
    exit(0);
  }

  if  ( int(2*J)%2 ==0){
    if ( int(2*lam1)%2 !=0 || int(2*lam1)%2 !=0){
      Alert << " J=\t" << J << "\tlam1=\t" << lam1 << "\tlam2=\t" << lam2 << "\t is not allowed" << endmsg;
    }
  }

  if  ( int(2*J)%2 !=0){
    if ( int(2*lam1)%2 ==0 || int(2*lam1)%2 ==0){
      Alert << " J=\t" << J << "\tlam1=\t" << lam1 << "\tlam2=\t" << lam2 << "\t is not allowed" << endmsg;
    }
  }

   std::stringstream Jstrstr;
   Jstrstr << _J;

   std::stringstream lam1strstr;
   lam1strstr << _lam1;

   std::stringstream lam2strstr;
   lam2strstr << _lam2;

//    std::string rootFileName="./Dfunc_J"+Jstrstr.str()+"_lam1_"+lam1strstr.str()+"_lam2_"+lam2strstr.str()+".root";
//    _theTFile=new TFile(rootFileName.c_str(),"recreate");

    _theTFile=new TFile("DfuncHist.root","recreate");
   std::string histName="d("+Jstrstr.str()+","+lam1strstr.str()+","+lam2strstr.str()+")";
  _dcostHist = new TH1F("_dcostHist",histName.c_str(),100, -1., 1.);
  _dcostHist->SetXTitle("cos(#Theta)"); 

  for (double cost=-1.; cost<1.; cost+=.02){
    double theta = acos(cost);
    complex<double> currentd = Wigner_D(0., theta,0, _J, _lam1, _lam2);
    // double weight=norm(currentd);
    double weight=currentd.real();
    _dcostHist->Fill( cost, weight );    
  }
}

DfuncHist::~DfuncHist()
{
   _theTFile->Write();
   _theTFile->Close();
}

