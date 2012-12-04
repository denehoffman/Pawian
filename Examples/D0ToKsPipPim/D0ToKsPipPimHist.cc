#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>

#include "Examples/D0ToKsPipPim/D0ToKsPipPimHist.hh"
#include "Examples/D0ToKsPipPim/D0ToKsPipPimEventList.hh"
#include "PwaUtils/KinUtils.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TLorentzVector.h"
#include "ErrLogger/ErrLogger.hh"


D0ToKsPipPimHist::D0ToKsPipPimHist(boost::shared_ptr<const EvtDataBaseListNew> theEvtList) :
  _theTFile(0),
  _dalitzDataHist(0),
  _dalitzMcHist(0),
  _dalitzFittedHist(0),
  _pippimMassDataHist(0),
  _pippimMassMcHist(0),
  _pippimMassFittedHist(0),
  _KspipMassDataHist(0),
  _KspipMassMcHist(0),
  _KspipMassFittedHist(0),
  _KspimMassDataHist(0),
  _KspimMassMcHist(0),
  _KspimMassFittedHist(0)
{
  if(0==theEvtList){
    Alert <<"D0ToKsPipPimEventList* theEvtList is a 0 pointer !!!!" ;  // << endmsg;
    exit(1);
  }

  initRootStuff();
  
  const std::vector<EvtDataNew*> dataList=theEvtList->getDataVecs();

  std::vector<EvtDataNew*>::const_iterator it=dataList.begin();
  while(it!=dataList.end())
    {
      const double evtWeight=(*it)->evtWeight;
      plotDalitz(_dalitzDataHist, (*it), evtWeight);
      plotPiPi(_pippimMassDataHist, (*it), evtWeight  );
      plotKsPip(_KspipMassDataHist, (*it), evtWeight );
      plotKsPim(_KspimMassDataHist, (*it), evtWeight );
      ++it;
    }

 const std::vector<EvtDataNew*> mcList=theEvtList->getMcVecs();  
  it=mcList.begin();
  while(it!=mcList.end())
    { 
     plotDalitz(_dalitzMcHist, (*it), 1.);
      plotPiPi(_pippimMassMcHist, (*it), 1.  );
      plotKsPip(_KspipMassMcHist, (*it), 1. );
      plotKsPim(_KspimMassMcHist, (*it), 1. );
      ++it;
    }

}

D0ToKsPipPimHist::D0ToKsPipPimHist(boost::shared_ptr<AbsLhNew> theLh, fitParamsNew& theFitParms):
  _theTFile(0),
  _dalitzDataHist(0),
  _dalitzMcHist(0),
  _dalitzFittedHist(0),
  _pippimMassDataHist(0),
  _pippimMassMcHist(0),
  _pippimMassFittedHist(0),
  _KspipMassDataHist(0),
  _KspipMassMcHist(0),
  _KspipMassFittedHist(0),
  _KspimMassDataHist(0),
  _KspimMassMcHist(0),
  _KspimMassFittedHist(0)
{ 
  if(0==theLh){
    Alert <<"AbsLh* is a 0 pointer !!!!" ;  // << endmsg;
    exit(1);
  }

   initRootStuff();
   boost::shared_ptr<const EvtDataBaseListNew> theEvtList=theLh->getEventList();
   const std::vector<EvtDataNew*> dataList=theEvtList->getDataVecs();  

   std::vector<EvtDataNew*>::const_iterator it=dataList.begin();
  while(it!=dataList.end())
    {
      const double evtWeight=(*it)->evtWeight;
      plotDalitz(_dalitzDataHist, (*it), evtWeight);
      plotPiPi(_pippimMassDataHist, (*it), evtWeight  );
      plotKsPip(_KspipMassDataHist, (*it), evtWeight );
      plotKsPim(_KspimMassDataHist, (*it), evtWeight );      
      ++it;
    }

 const std::vector<EvtDataNew*> mcList=theEvtList->getMcVecs();  
  it=mcList.begin();
  while(it!=mcList.end())
    { 
      plotDalitz(_dalitzMcHist, (*it), 1.);
      plotPiPi(_pippimMassMcHist, (*it), 1.  );
      plotKsPip(_KspipMassMcHist, (*it), 1. );
      plotKsPim(_KspimMassMcHist, (*it), 1. );      

      double evtWeight= theLh->calcEvtIntensity( (*it), theFitParms );
      plotDalitz(_dalitzFittedHist, (*it), evtWeight);
      plotPiPi(_pippimMassFittedHist, (*it), evtWeight  );
      plotKsPip(_KspipMassFittedHist, (*it), evtWeight );
      plotKsPim(_KspimMassFittedHist, (*it), evtWeight ); 
       ++it;
     }

   double integralDataWoWeight=(double) theEvtList->getDataVecs().size();
   Info <<"No of data events without weight " << integralDataWoWeight ;  // << endmsg;

   double integralDataWWeight=(double) theEvtList->NoOfWeightedDataEvts();
   Info <<"No of data events with weight " << integralDataWWeight ;  // << endmsg;


   double integralMC=(double) theEvtList->NoOfWeightedMcEvts();
   Info <<"No of MC events " << integralMC ;  // << endmsg; 

   Info <<"scaling factor  " << integralDataWWeight/integralMC ;  // << endmsg;

   double scaleFactor = integralDataWWeight/integralMC;
   _dalitzFittedHist->Scale(scaleFactor);
   _pippimMassFittedHist->Scale(scaleFactor);
   _KspipMassFittedHist->Scale(scaleFactor);
   _KspimMassFittedHist->Scale(scaleFactor);

   double integralFitted=(double) _pippimMassFittedHist->Integral();
   Info <<"No of fit events " << integralFitted ;  // << endmsg; 
 }


D0ToKsPipPimHist::~D0ToKsPipPimHist()
{
  _theTFile->Write();
  _theTFile->Close();
}

void D0ToKsPipPimHist::initRootStuff()
{
  std::string rootFileName="./D0ToKsPipPim.root";

  _theTFile=new TFile(rootFileName.c_str(),"recreate");
  int xbins = 200;
  double xmin=0.;
  double xmax =3.2;
  int ybins=xbins;
  double ymin=xmin;
  double ymax=xmax;
  
  _dalitzDataHist= new TH2F("_dalitzDataHist","Dpl K_{s} #pi^{+} #pi^{-} data",xbins, xmin, xmax, ybins, ymin, ymax );
  _dalitzMcHist= new TH2F("_dalitzMcHist","Dpl K+K- K+#pi^{0} MC",xbins, xmin, xmax, ybins, ymin, ymax);
  _dalitzFittedHist= new TH2F("_dalitzFittedHist","Dpl K+K- K+#pi^{0} fit",xbins, xmin, xmax, ybins, ymin, ymax );

  
  int nbins=200;
   xmin=0.2;
   xmax=1.5;
  _pippimMassDataHist = new TH1F("_pippimMassDataHist", "M(#pi^{+} #pi^{-}) data", nbins, xmin, xmax);
  _pippimMassMcHist = new TH1F("_pippimMassMcHist", "M(#pi^{+} #pi^{-}) MC", nbins, xmin, xmax);
  _pippimMassFittedHist = new TH1F("_pippimMassFittedHist", "M(#pi^{+} #pi^{-}) fit", nbins, xmin, xmax);

  nbins=200;  
  xmin=0.5;
  xmax=1.8;
  _KspipMassDataHist = new TH1F("_KspipMassDataHist", "M(K_{s} #pi^{+}) data", nbins, xmin, xmax);
  _KspipMassMcHist = new TH1F("_KspipMassMcHist", "M(K_{s} #pi^{+}) MC", nbins, xmin, xmax);
  _KspipMassFittedHist = new TH1F("_KspipMassFittedHist", "M(K_{s} #pi^{+}) fit", nbins, xmin, xmax);

  _KspimMassDataHist = new TH1F("_KspimMassDataHist", "M(K_{s} #pi^{-}) data", nbins, xmin, xmax);
  _KspimMassMcHist = new TH1F("_KspimMassMcHist", "M(K_{s} #pi^{-}) MC", nbins, xmin, xmax);
  _KspimMassFittedHist = new TH1F("_KspimMassFittedHist", "M(K_{s} #pi^{-}) fit", nbins, xmin, xmax);

}

void D0ToKsPipPimHist::plotDalitz(TH2F* theHisto,  EvtDataNew* theData, double weight)
{
  Vector4<double>& V4_KsPip_HeliD0 = theData->FourVecsDec[enumD0KsPiPi4V::V4_KsPip_HeliD0]  ;
  Vector4<double>& V4_KsPim_HeliD0 = theData->FourVecsDec[enumD0KsPiPi4V::V4_KsPim_HeliD0]  ;
  theHisto->Fill( V4_KsPip_HeliD0.M()*V4_KsPip_HeliD0.M(), V4_KsPim_HeliD0.M()*V4_KsPim_HeliD0.M(),weight);  
}

void D0ToKsPipPimHist::plotPiPi(TH1F* theHisto, EvtDataNew* theData, double weight){
  Vector4<double> v4 = theData->FourVecsDec[enumD0KsPiPi4V::V4_PipPim_HeliD0];
  theHisto->Fill( v4.M(), weight );
}

void D0ToKsPipPimHist::plotKsPip(TH1F* theHisto, EvtDataNew* theData, double weight){
  Vector4<double>& v4 = theData->FourVecsDec[enumD0KsPiPi4V::V4_KsPip_HeliD0];
  theHisto->Fill( v4.M(), weight ); 
}

void D0ToKsPipPimHist::plotKsPim(TH1F* theHisto, EvtDataNew* theData, double weight){
  Vector4<double>& v4 = theData->FourVecsDec[enumD0KsPiPi4V::V4_KsPim_HeliD0];
  theHisto->Fill( v4.M(), weight ); 
}





