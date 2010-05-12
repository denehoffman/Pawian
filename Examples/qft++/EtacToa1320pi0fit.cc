#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/qft++/EtacToa1320pi0fit.hh"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"


EtacToa1320pi0fit::EtacToa1320pi0fit(int kindOfData) :
  _kindOfData(kindOfData),
  _nOfData(10000),
  _nOfMc(30000),
  _fitHistsFilled(false),
  _theTFile(0),
  _dalitzDataHist(0),
  _dalitzMcHist(0),
  _dalitzFittedHist(0),
  _invpietaDataHist(0),
  _invpietaMcHist(0),
  _invpietaFittedHist(0),
  _invpipiDataHist(0),
  _invpipiMcHist(0),
  _invpipiFittedHist(0)
{
  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR");
  _mcPath=theSourcePath+std::string("/Examples/qft++/data/mc100Mgev.dat");

  if (kindOfData==0) _dataPath=theSourcePath+std::string("/Examples/qft++/data/dataEtacToA0Pi_100Mgev.dat");
  else if (kindOfData==2)   _dataPath=theSourcePath+std::string("/Examples/qft++/data/dataSpin2100MgevNew.dat");
  else{
    std::cout <<"this kind of data: " << kindOfData << "  is not supported!!!!\n"
	      <<"initialze EtacToa1320pi0fit either with 0 or with 2 !!!!"
	      << std::endl;
    assert(0);    
  }

  initRootStuff();
  read4Vecs(_dataPath, _nOfData, _data4Vecs);
  read4Vecs(_mcPath,_nOfMc, _mc4Vecs);

  for (size_t it=0; it<_data4Vecs.size(); it++)
    { 
      plotDalitz(_dalitzDataHist, _data4Vecs[it], 1);
      plotInvPiEta(_invpietaDataHist, _data4Vecs[it], 1);
      plotInvPiPi(_invpipiDataHist, _data4Vecs[it], 1);
    }

  for (size_t it=0; it<_mc4Vecs.size(); it++)
    { 
      plotDalitz(_dalitzMcHist, _mc4Vecs[it], 1);
      plotInvPiEta(_invpietaMcHist, _mc4Vecs[it], 1);
      plotInvPiPi(_invpipiMcHist, _mc4Vecs[it], 1);
    } 


}



EtacToa1320pi0fit::~EtacToa1320pi0fit()
{
  _theTFile->Write();
  _theTFile->Close();
//   if (0!=_theTFile) delete _theTFile;
//   if (0!=_dalitzDataHist) delete _dalitzDataHist;
//   if (0!=_dalitzMcHist) delete _dalitzMcHist;
//   if (0!=_dalitzFittedHist) delete _dalitzFittedHist;
//   if (0!=_invpietaDataHist) delete _invpietaDataHist;
//   if (0!=_invpietaMcHist) delete _invpietaMcHist;
//   if (0!=_invpietaFittedHist) delete _invpietaFittedHist;
//   if (0!=_invpipiDataHist) delete _invpipiDataHist;
//   if (0!=_invpipiMcHist) delete _invpipiMcHist;
//   if (0!=_invpipiFittedHist) delete _invpipiFittedHist;
}

void EtacToa1320pi0fit::initRootStuff()
{
  _theTFile=new TFile("./a1320Fit.root","recreate");
  _dalitzDataHist= new TH2F("_dalitzDataHist","dalitzplot data",50, 0., 8.,50, 0., 8.);
  _dalitzMcHist= new TH2F("_dalitzMcHist","dalitzplot mc",50, 0., 8.,50, 0., 8.);
  _dalitzFittedHist= new TH2F("_dalitzFittedHist","dalitzplot fit",50, 0., 8.,50, 0., 8.);
  _invpietaDataHist= new TH1F("_invpietaDataHist","M(#pi #eta) data",120, 0.6, 3.0);
  _invpietaMcHist= new TH1F("_invpietaMcHist","M(#pi #eta) mc",120, 0.6, 3.0);
  _invpietaFittedHist= new TH1F("_invpietaFittedHist","M(#pi #eta) fit",120, 0.6, 3.0);
  _invpipiDataHist= new TH1F("_invpipiDataHist","M(#pi #pi) data",120, 0.2, 2.6);
  _invpipiMcHist= new TH1F("_invpipiMcHist","M(#pi #pi) mc",120, 0.2, 2.6);
  _invpipiFittedHist= new TH1F("_invpipiFittedHist","M(#pi #pi) fit",120, 0.2, 2.6);
}



void EtacToa1320pi0fit::read4Vecs(std::string& path,  int nEvts, std::vector<evt4Vec>& the4Vecs)
{
  std::cout << "calculate 4Vecs and amplitudes for " << path << std::endl;
  std::ifstream inputStream(path.c_str(), std::ios::in);

   if (!inputStream.good()) 
   {
     std::cout << "Input file doesn't exsits !!! ";
     assert(0);
   }

   int counter=0; 
   while (!inputStream.eof() && counter<nEvts )
     {
       if ( counter%1000 == 0 ) std::cout << "event " << counter << std::endl;
       Vector4<double> pi14V,pi24V,eta4V; // 4-momenta
       get4Vecs(inputStream, pi14V);
       get4Vecs(inputStream, pi24V);
       get4Vecs(inputStream, eta4V);

       Vector4<double>  cm_4V(pi14V+pi24V+eta4V);

       Vector4<double> pi1_cm_4V(pi14V);
       pi1_cm_4V.Boost(cm_4V);
       Vector4<double>  pi2_cm_4V(pi24V);
       pi2_cm_4V.Boost(cm_4V);
       Vector4<double>  eta_cm_4V(eta4V);
       eta_cm_4V.Boost(cm_4V);

       Vector4<double> a2_0_cm_4Vec(pi1_cm_4V+eta_cm_4V);
       Vector4<double> a2_1_cm_4Vec(pi2_cm_4V+eta_cm_4V);
       Vector4<double>  cm_cm_4V(cm_4V);
       cm_cm_4V.Boost(cm_4V);

       evt4Vec theEvt4Vec;
       theEvt4Vec.pi0_0_4Vec=pi1_cm_4V;
       theEvt4Vec.pi0_1_4Vec=pi2_cm_4V;
       theEvt4Vec.eta_4Vec=eta_cm_4V;
       theEvt4Vec.a2_0_4Vec=a2_0_cm_4Vec;
       theEvt4Vec.a2_1_4Vec=a2_1_cm_4Vec;
       theEvt4Vec.cm_4Vec=cm_cm_4V;

       theEvt4Vec.spin2_0_amp = calcSpin2Amp(theEvt4Vec.pi0_0_4Vec, theEvt4Vec.eta_4Vec, theEvt4Vec.pi0_1_4Vec);

       theEvt4Vec.spin2_1_amp = calcSpin2Amp(theEvt4Vec.pi0_1_4Vec, theEvt4Vec.eta_4Vec, theEvt4Vec.pi0_0_4Vec);

       theEvt4Vec.spin1_0_amp = calcSpin1Amp(theEvt4Vec.pi0_0_4Vec, theEvt4Vec.eta_4Vec, theEvt4Vec.pi0_1_4Vec);

       theEvt4Vec.spin1_1_amp = calcSpin1Amp(theEvt4Vec.pi0_1_4Vec, theEvt4Vec.eta_4Vec, theEvt4Vec.pi0_0_4Vec);

       the4Vecs.push_back(theEvt4Vec);
       counter++;
     } 
}


void EtacToa1320pi0fit::get4Vecs (std::ifstream& inStream, Vector4<double>& the4Vec)
{
  double tmpPx,tmpPy, tmpPz, tmpE;
  inStream >> tmpPx;
  inStream >> tmpPy;
  inStream >> tmpPz;
  inStream  >> tmpE;
  the4Vec.SetP4(tmpE, tmpPx, tmpPy, tmpPz);  
}


void EtacToa1320pi0fit::plotDalitz(TH2F* theHisto, evt4Vec& theEvt4Vecs, double weight)
{

      Vector4<double> pi0eta4V=theEvt4Vecs.pi0_0_4Vec + theEvt4Vecs.eta_4Vec;
      Vector4<double> pi1eta4V=theEvt4Vecs.pi0_1_4Vec + theEvt4Vecs.eta_4Vec;
      theHisto->Fill(pi0eta4V.M()*pi0eta4V.M(),pi1eta4V.M()*pi1eta4V.M(), weight);
      theHisto->Fill(pi1eta4V.M()*pi1eta4V.M(),pi0eta4V.M()*pi0eta4V.M(), weight);

}

void EtacToa1320pi0fit::plotInvPiEta (TH1F* theHisto, evt4Vec& theEvt4Vecs, double weight)
{
      Vector4<double> pi0eta4V=theEvt4Vecs.pi0_0_4Vec + theEvt4Vecs.eta_4Vec;
      Vector4<double> pi1eta4V=theEvt4Vecs.pi0_1_4Vec + theEvt4Vecs.eta_4Vec;
      theHisto->Fill(pi0eta4V.M(), weight);
      theHisto->Fill(pi1eta4V.M(), weight);

}

void EtacToa1320pi0fit::plotInvPiPi (TH1F* theHisto, evt4Vec& theEvt4Vecs, double weight)
{

      Vector4<double> pi0pi0=theEvt4Vecs.pi0_0_4Vec + theEvt4Vecs.pi0_1_4Vec;
      theHisto->Fill(pi0pi0.M(), weight);

}


Tensor<complex<double> > EtacToa1320pi0fit::calcSpin2Amp(Vector4<double>& pi_cm, Vector4<double>& eta_cm, Vector4<double>& pi_recoil_cm)
{

  Tensor<complex<double> > result(0); //rank-0 Tensor for the amplitude

  Vector4<double> inter_cm=eta_cm + pi_cm;
  
  PolVector inter_Pol_cm(2);
  inter_Pol_cm.SetP4(inter_cm,inter_cm.M());

  OrbitalTensor orb_eta_To_inter_piRecoil(2);
  orb_eta_To_inter_piRecoil.SetP4(inter_cm, pi_recoil_cm);
  
  OrbitalTensor orb_inter_To_pieta(2);
  orb_inter_To_pieta.SetP4(pi_cm, eta_cm);

  result = (orb_inter_To_pieta) | ( inter_Pol_cm.Projector() | orb_eta_To_inter_piRecoil );
 
  return result;
}

Tensor<complex<double> > EtacToa1320pi0fit::calcSpin1Amp(Vector4<double>& pi_cm, Vector4<double>& eta_cm, Vector4<double>& pi_recoil_cm)
{


  Tensor<complex<double> > result(0); // rank-0 Tensor for the amplitude

  Vector4<double> inter_cm=eta_cm + pi_cm;
  
  PolVector inter_Pol_cm(1);
  inter_Pol_cm.SetP4(inter_cm,inter_cm.M());

  OrbitalTensor orb_eta_To_inter_piRecoil(1);
  orb_eta_To_inter_piRecoil.SetP4(inter_cm, pi_recoil_cm);
  
  OrbitalTensor orb_inter_To_pieta(1);
  orb_inter_To_pieta.SetP4(pi_cm, eta_cm);

  result = (orb_inter_To_pieta) | ( inter_Pol_cm.Projector() | orb_eta_To_inter_piRecoil );
 
  return result;
}

double EtacToa1320pi0fit::calcIntensityCache(evt4Vec& theEvtVec, const fitParamVal& theParamVal)
{
  double spin2cont=theParamVal.cont2spin;
  double spin1cont=theParamVal.cont1spin;
  double spin0cont=theParamVal.cont0spin;
  
  double intensity=0;

  Tensor<complex<double> > amp(0); // rank-0 Tensor for the amplitude

  amp = theEvtVec.spin2_0_amp*BreitWigner(theEvtVec.a2_0_4Vec, theParamVal.a1320Mass, theParamVal.a1320Width);
  intensity += spin2cont*norm(amp(0,0));

  amp = theEvtVec.spin2_1_amp*BreitWigner(theEvtVec.a2_1_4Vec, theParamVal.a1320Mass, theParamVal.a1320Width);
  intensity += spin2cont*norm(amp(0,0));

  amp = theEvtVec.spin1_0_amp*BreitWigner(theEvtVec.a2_0_4Vec, theParamVal.a1320Mass, theParamVal.a1320Width);
  intensity += spin1cont*norm(amp(0,0));

  amp = theEvtVec.spin1_1_amp*BreitWigner(theEvtVec.a2_1_4Vec, theParamVal.a1320Mass, theParamVal.a1320Width);
  intensity += spin1cont*norm(amp(0,0));

  //dont know how to normalize
  amp = 0.1*BreitWigner(theEvtVec.a2_0_4Vec, theParamVal.a1320Mass, theParamVal.a1320Width);
  intensity += spin0cont*norm(amp(0,0));

  amp = 0.1*BreitWigner(theEvtVec.a2_1_4Vec, theParamVal.a1320Mass, theParamVal.a1320Width);
  intensity += spin0cont*norm(amp(0,0));


  return intensity;
}



double EtacToa1320pi0fit::calcLogLh(const fitParamVal& theParamVal)
{
  double logLH=0.;
  double logLH_data=0.;
  for (size_t it=0; it<_data4Vecs.size(); it++)
    { 
      double intensity=calcIntensityCache(_data4Vecs[it], theParamVal);
      if (intensity>0.) logLH_data+=log10(intensity);
    }

  std::cout << "logLH_data= " << logLH_data << std::endl;


  double LH_mc=0.;
  for (size_t it=0; it<_mc4Vecs.size(); it++)
    { 
      double intensity=calcIntensityCache(_mc4Vecs[it], theParamVal);
      LH_mc+=intensity;
    }
  std::cout << "LH_mc= " << LH_mc << std::endl;

  double logLH_mc_Norm=0.;
  if (LH_mc>0.) logLH_mc_Norm=log10(LH_mc/_mc4Vecs.size());


  logLH=_data4Vecs.size()/2.*(LH_mc/_mc4Vecs.size()-1)*(LH_mc/_mc4Vecs.size()-1)
    -logLH_data
    +_data4Vecs.size()*logLH_mc_Norm;

 return logLH;
}

bool EtacToa1320pi0fit::initFitParameters(MnUserParameters& minuitParams)
{
  bool result=true;
  if (_kindOfData==2)
    {  
      minuitParams.Add("InterMass", 1.6, .1, 2.3, 0.6);
      minuitParams.Add("InterWidth", 0.04, .01, 0.8, 0.01);
      minuitParams.Add("spin0", 0.3, .1, 1., 0.);
      minuitParams.Add("spin1", 0.3, .1, 1., 0.);
      minuitParams.Add("spin2", 0.3, .1, 1., 0.);
    }
  else  //_kindOfData==0
    {
      minuitParams.Add("InterMass", 1.1, .1, 2.3, 0.6);
      minuitParams.Add("InterWidth", 0.04, .01, 0.8, 0.01);
      minuitParams.Add("spin0", 0.3, .1, 1., 0.);
      minuitParams.Add("spin1", 0.3, .1, 1., 0.);
      minuitParams.Add("spin2", 0.3, .1, 1., 0.);
    }
  
  return result;
}

bool EtacToa1320pi0fit::setFitParamVal(fitParamVal& fitParamVal, const std::vector<double>& par)
{
  bool result=true;
  fitParamVal.a1320Mass=par[0];
  fitParamVal.a1320Width=par[1];
  fitParamVal.cont0spin=par[2];
  fitParamVal.cont1spin=par[3];
  fitParamVal.cont2spin=par[4];
  return result;
} 

bool  EtacToa1320pi0fit::fillFitHists(const fitParamVal& fitParamVal)
{
  if(_fitHistsFilled)
    {
      std::cout <<"Not possible to fill the histos with fit results!!!!\n"
		<<"They are already filled!!!!"
		<<std::endl;
      return false;
    }

  for (size_t it=0; it<_mc4Vecs.size(); it++)
    { 
      double evtWeight=calcIntensityCache(_mc4Vecs[it], fitParamVal);
      plotDalitz(_dalitzFittedHist, _mc4Vecs[it], evtWeight);
      plotInvPiEta(_invpietaFittedHist, _mc4Vecs[it], evtWeight);
      plotInvPiPi(_invpipiFittedHist, _mc4Vecs[it], evtWeight);
    }
  //normalize fitted histos

  double integralData=_invpipiDataHist->Integral();
  std::cout << "integralData= " << integralData << std::endl;
  double integralFittedMc=_invpipiFittedHist->Integral();  
  std::cout << "integralFittedMc= " << integralFittedMc << std::endl;
  _invpipiFittedHist->Scale(integralData/integralFittedMc);
  _invpietaFittedHist->Scale(integralData/integralFittedMc);
  _dalitzFittedHist->Scale(integralData/integralFittedMc);

  return true;

}

