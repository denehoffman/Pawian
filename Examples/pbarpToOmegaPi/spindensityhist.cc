#include "Examples/pbarpToOmegaPi/spindensityhist.hh"
#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
//#include "Examples/pbarpToOmegaPi/OmegaPiHist.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiEventList.hh"
//#include "Examples/pbarpToOmegaPi/OmegaPiLh.hh"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "ErrLogger/ErrLogger.hh"
#include "Examples/pbarpToOmegaPi/spindensity.hh"
#include <cstdlib>

SpinDensityHist::SpinDensityHist(const std::string &thePathToRootFile, 
                                 const std::vector<OmegaPiData::OmPiEvtData*> &theEventData,
                                 OmegaPiData::fitParamVal &theParamVal):
  _numOfEvts(1000),
  m_EventData(theEventData),
  m_PathToRootFile(thePathToRootFile)
{
  m_pfitParamVal = &theParamVal;
  _theTFile=new TFile(m_PathToRootFile.c_str(),"recreate");
  
  _cosOmegaHeliMcHist= new TH1F("_cosOmegaHeliMcHist","cos(#Theta) #omega heli MC",101, -1.0, 1.0);

  unsigned int counter=0;
  std::vector<OmPiEvtData*>::const_iterator it=m_EventData.begin();
  while(it!=m_EventData.end())
  {
    if (counter>_numOfEvts) break;
     _cosOmegaHeliMcHist->Fill((*it)->omegaHeliCm4Vec.CosTheta(), 1.0);
    ++it;
    counter++;
  }
}


SpinDensityHist::~SpinDensityHist()
{
  _theTFile->Write();
  _theTFile->Close();
}

void SpinDensityHist::createHistogram(int M, int M_)
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
  createSpinDensityHist(_spinDensityDataHist,  M, M_);
  
  //Create spin density histogram with imaginary part
  strstreamName.str("");
  strstreamName << "_" << "Imag" << strCohIncoh << "spinDensityDataHist_M1" << strPrefixM1 << Mn << "M2" << strPrefixM2 << Mn_;
  TH1F* _spinDensityDataHistImag= new TH1F(strstreamName.str().c_str(),strstreamTitle.str().c_str(),101, -1.0, 1.0);
  createSpinDensityHist(_spinDensityDataHistImag,  M, M_);

  _spinDensityDataHist->Divide(_cosOmegaHeliMcHist);
  _spinDensityDataHistImag->Divide(_cosOmegaHeliMcHist);

  DebugMsg << "calculation done for spin density matrix element rho" << M << M_ << endmsg;
}

void SpinDensityHist::createHistograms()
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

void SpinDensityHist::AddData(TH1F* theHisto, const OmPiEvtData& theEvtData, double dSpinDensity)
{
  if (theHisto != 0)  
  {
       theHisto->Fill(theEvtData.omegaHeliCm4Vec.CosTheta(), dSpinDensity);
  }  
}


/*!
    \fn SpinDensityHist::createSpinDensityHist()
 */
void SpinDensityHist::createSpinDensityHist(TH1F* theHisto, int M, int M_,bool bReal)
{
        
  std::vector<OmegaPiData::OmPiEvtData*>::const_iterator iterd;
  complex<double> SpinDensity;
  unsigned int counter=0;

  for (iterd=m_EventData.begin(); iterd!=m_EventData.end(); ++iterd)
  {
    if (counter>_numOfEvts) continue;
    SpinDensity=spinDensity::calcSpinDensityIncoherent(M, M_ , (*iterd), (*m_pfitParamVal));
    //    else SpinDensity=spinDensity::calcSpinDensityCoherent(M, M_ , (*iterd), (*m_pfitParamVal));

//     if (bIncoherent) SpinDensity=spinDensity::calcSpinDensityOmegaFrame(M, M_ , (*iterd), (*m_pfitParamVal));
//     else SpinDensity=spinDensity::calcSpinDensityCoherent(M, M_ , (*iterd), (*m_pfitParamVal)); 
    
    
    if (bReal) AddData(theHisto, *(*iterd), SpinDensity.real());
    else AddData(theHisto, *(*iterd), SpinDensity.imag());
    counter++;
  } 
        
}
