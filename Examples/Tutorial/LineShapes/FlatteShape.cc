#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Examples/Tutorial/LineShapes/FlatteShape.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"


#include "ErrLogger/ErrLogger.hh"

FlatteShape::FlatteShape(std::string ptype, double g1, double g2) :
  _theTFile(0)
  ,_histShapeLow(0)
  ,_histShapeHigh(0)
  ,_pType(ptype)
{
  if(_pType!="a0" && _pType!="f0"){
    Alert <<"only a0 and f0 are supported" << endmsg;
    exit(1);
  }
  if(g1<0.  || g2<0.){
    Alert <<"requirement: g1>0. and g2>0." << endmsg;
    exit(1);
  }

  std::string rootFileName="./FlatteShape"+_pType+".root";
  _theTFile=new TFile(rootFileName.c_str(),"recreate");


  double mass1=0.;
  double mass2=0.;

  if(_pType=="a0"){
    mass1 = 0.1349766;
    mass2 = 0.547853;
  }
  else if(_pType=="f0"){
    mass1 = 0.1349766;
    mass2= 0.1349766;
  }

  std::pair <const double, const double> decPairLow=make_pair(mass1, mass2); 
  const double KplusMass = 0.493677;
  const double K0Mass = 0.497614;
  std::pair <const double, const double> decPairHigh=make_pair(KplusMass,K0Mass);


  int size=700;
  double massMin=0.35;
  double massMax=1.85;

//   double massMin=1.;
//   double massMax=1.004;
  double stepSize=(massMax-massMin)/size;

  _histShapeLow= new TH1F("_histShapeLow","hist low",size+1, massMin, massMax);
  _histShapeHigh= new TH1F("_histShapeHigh","hist high",size+1, massMin, massMax);

  for (double mass=massMin; mass<massMax; mass+=stepSize){
    Vector4<double> mass4Vec(mass, 0.,0.,0.);

    _histShapeLow->Fill(mass4Vec.M(), norm(Flatte(mass4Vec, decPairLow, decPairHigh, 0.981, g1, g2)));
    _histShapeHigh->Fill(mass4Vec.M(), norm(Flatte(mass4Vec, decPairHigh, decPairLow, 0.981, g2, g1))); 
  }

}

FlatteShape::~FlatteShape()
{
  _theTFile->Write();
  _theTFile->Close();
}

