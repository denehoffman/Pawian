#include "Examples/pbarpToOmegaPiLS/SpinDensityHistLS.hh"
#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Examples/pbarpToOmegaPiLS/AbsOmegaPiLhLS.hh"
#include "Examples/pbarpToOmegaPiLS/AbsOmegaPiEventListLS.hh"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "ErrLogger/ErrLogger.hh"
#include <cstdlib>

SpinDensityHistLS::SpinDensityHistLS(const std::string &thePathToRootFile, 
				     boost::shared_ptr<AbsOmegaPiLhLS> absOmegaPiLh,
				     OmegaPiDataLS::fitParamVal &theParamVal):
  _numOfEvts(1000),
  _omegaPiLh(absOmegaPiLh),
  m_PathToRootFile(thePathToRootFile)
{
  boost::shared_ptr<const AbsOmegaPiEventListLS> theEvtList=_omegaPiLh->getEventList();
  m_EventData=theEvtList->getMcVecs();
  m_pfitParamVal = &theParamVal;
  _theTFile=new TFile(m_PathToRootFile.c_str(),"recreate");
  
  _cosOmegaHeliMcHist= new TH1F("_cosOmegaHeliMcHist","cos(#Theta) #omega heli MC",101, -1.0, 1.0);

  unsigned int counter=0;
  std::vector<OmPiEvtDataLS*>::const_iterator it=m_EventData.begin();
  while(it!=m_EventData.end())
  {
    if (counter>_numOfEvts) break;
     _cosOmegaHeliMcHist->Fill((*it)->omegaHeliCm4Vec.CosTheta(), 1.0);
    ++it;
    counter++;
  }
}


SpinDensityHistLS::~SpinDensityHistLS()
{
  _theTFile->Write();
  _theTFile->Close();
}

void SpinDensityHistLS::createHistogram(int M, int M_)
{

  std::stringstream strstreamTitle;
  strstreamTitle << "#rho^{0}_{" << M << " " << M_ << "}";
  std::stringstream strstreamName;
  string strPrefixM1;
  string strPrefixM2;
  int Mn=0, Mn_=0;
  if (M<0) 
  {
    strPrefixM1="m";
    Mn=abs(M);
  }
  else Mn=M;
  
  if (M_<0) 
  {
    strPrefixM2="m";
    Mn_=abs(M_);
  }
  else Mn_=M_;
  
  string strCohIncoh;
  strCohIncoh="Incoh";

  //Create spin density histogram with real part
  strstreamName << "_" << strCohIncoh << "spinDensityDataHist_M1" << strPrefixM1 << Mn << "M2" << strPrefixM2 << Mn_;
  TH1F* _spinDensityDataHist= new TH1F(strstreamName.str().c_str(),strstreamTitle.str().c_str(),101, -1.0, 1.0);
  createSpinDensityHist(_spinDensityDataHist,  M, M_, true);
  
  //Create spin density histogram with imaginary part
  strstreamName.str("");
  strstreamName << "_" << "Imag" << strCohIncoh << "spinDensityDataHist_M1" << strPrefixM1 << Mn << "M2" << strPrefixM2 << Mn_;
  TH1F* _spinDensityDataHistImag= new TH1F(strstreamName.str().c_str(),strstreamTitle.str().c_str(),101, -1.0, 1.0);
  createSpinDensityHist(_spinDensityDataHistImag,  M, M_, false);

  _spinDensityDataHist->Divide(_cosOmegaHeliMcHist);
  _spinDensityDataHistImag->Divide(_cosOmegaHeliMcHist);

  DebugMsg << "calculation done for spin density matrix element rho" << M << M_ << endmsg;
}

void SpinDensityHistLS::createHistograms()
{
  createHistogram(0,0);
  createHistogram(1,0);
  createHistogram(0,1);
  createHistogram(-1,0);
  createHistogram(0,-1);
  createHistogram(1,1);
  createHistogram(1,-1);
  createHistogram(-1,1);
  createHistogram(-1,-1);
  
}

void SpinDensityHistLS::AddData(TH1F* theHisto, const OmPiEvtDataLS& theEvtData, double dSpinDensity)
{
  if (theHisto != 0)  
  {
       theHisto->Fill(theEvtData.omegaHeliCm4Vec.CosTheta(), dSpinDensity);
  }  
}


/*!
    \fn SpinDensityHist::createSpinDensityHist()
 */
void SpinDensityHistLS::createSpinDensityHist(TH1F* theHisto, int M, int M_,bool bReal)
{
        
  std::vector<OmegaPiDataLS::OmPiEvtDataLS*>::const_iterator iterd;
  complex<double> SpinDensity;
  unsigned int counter=0;

  boost::shared_ptr<const AbsOmegaPiEventListLS> theEvtList=_omegaPiLh->getEventList();
  std::vector<OmegaPiDataLS::OmPiEvtDataLS*> theData=theEvtList->getMcVecs();


//   for (iterd=m_EventData.begin(); iterd!=m_EventData.end(); ++iterd)
   for (iterd=theData.begin(); iterd!=theData.end(); ++iterd)
  {
    if (counter>_numOfEvts) continue;
    SpinDensity=_omegaPiLh->spinDensity(M, M_ , (*iterd), (*m_pfitParamVal));
//     SpinDensity=_omegaPiLh->spinDensityOmegaFrame(M, M_ , (*iterd), (*m_pfitParamVal));
    
    
    if (bReal) AddData(theHisto, *(*iterd), SpinDensity.real());
    else AddData(theHisto, *(*iterd), SpinDensity.imag());
    counter++;
  } 
        
}
