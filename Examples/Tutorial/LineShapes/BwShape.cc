#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Examples/Tutorial/LineShapes/BwShape.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"


#include "ErrLogger/ErrLogger.hh"

BwShape::BwShape(double MassRes, double MassWidth, double MassDec1, double MassDec2, unsigned int Lmax, double deltaMass) :
  _theTFile(0)
{

  if (MassRes <= MassDec1+MassDec2){
    Alert << "Mass of the resonance smaller than the masses of the decay particles\n" 
	  << MassRes << " < " << MassDec1+MassDec2 << endmsg;
    exit(0); 
  }

   std::stringstream Lmaxstrstr;
   Lmaxstrstr << Lmax;
   std::string rootFileName="./BwShapeLmax"+Lmaxstrstr.str()+".root";
   _theTFile=new TFile(rootFileName.c_str(),"recreate");
 
  for (unsigned int i=0; i<=Lmax; ++i){

    std::stringstream Lstrstr;
    Lstrstr << i;
    std::string histName="BreitWigner_L"+Lstrstr.str();


   _histMap[i]= new TH1F(histName.c_str(),histName.c_str(),301, MassRes-deltaMass, MassRes+deltaMass);

    histName="Argand_L"+Lstrstr.str();
    TH2F* currentArgandHist=new TH2F(histName.c_str(),histName.c_str(),301, -1., 1., 301, 0., 1.3);
    currentArgandHist->SetXTitle("Re(Bw)");
    currentArgandHist->SetYTitle("Im(Bw)");
    _argandHistMap[i]=currentArgandHist;
    
    
    //initialize here the phase shift histogramms _phaseHistMap
  }
  
  double stepSize=2.*deltaMass/300;
  
  for (unsigned int lIt=0; lIt<=Lmax; ++lIt){
    
    TH1F* currentHist=_histMap[lIt];
    TH2F* currentArgandHist=_argandHistMap[lIt];

    for (double massIt=MassRes-deltaMass; massIt<MassRes+deltaMass; massIt+=stepSize){
      
      Vector4<double>  res4V(massIt, 0., 0., 0.);
      
      complex<double> currentBW=BreitWignerBlattW(res4V, MassDec1, MassDec2, MassRes, MassWidth, lIt);
      double weight=norm(currentBW);
      currentHist->Fill(massIt,weight);
      
      currentArgandHist->Fill(currentBW.real(),currentBW.imag());
      
      //fill here the phase shift histogramms: hint "atan2(imag,real)"

    }

  }
}

BwShape::~BwShape()
{
   _theTFile->Write();
   _theTFile->Close();
}

