// pbarpEventList class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#include <sstream>

#include "pbarpUtils/spinDensityHist.hh"
#include "pbarpUtils/pbarpEnv.hh"
#include "pbarpUtils/pbarpBaseLh.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"

#include "TH1F.h"
#include "TFile.h"



spinDensityHist::spinDensityHist(boost::shared_ptr<AbsLh> theLh, fitParams& theFitParams) :
   _nBins(101)
  ,_maxEvents(2000)
  ,_theLh(theLh)
{
   _dataList=_theLh->getEventList()->getMcVecs();
   boost::dynamic_pointer_cast<pbarpBaseLh>(_theLh)->updateFitParams(theFitParams);

   std::stringstream spinDensityRootFileName;
   spinDensityRootFileName << "./spinDensity" << pbarpEnv::instance()->outputFileNameSuffix() << ".root";

   _spinDensityRootFile = new TFile(spinDensityRootFileName.str().c_str(), "recreate");

   std::vector<string> spinDensityParticles = pbarpEnv::instance()->spinDensityNames();
   std::vector<string>::iterator it;

   for(it=spinDensityParticles.begin(); it!=spinDensityParticles.end(); ++it){
      Info << "Calculating spin density matrix for particle " << (*it) << endmsg;
      calcSpinDensityMatrix(*it);
   }
}



spinDensityHist::~spinDensityHist(){
   _spinDensityRootFile->Write();
   _spinDensityRootFile->Close();
}



void spinDensityHist::calcSpinDensityMatrix(std::string& particleName){

   int J = pbarpEnv::instance()->particleTable()->particle(particleName)->J();

   if(J!=1){
      Alert << "Spin density calculation for J!=1 currently not supported." << endmsg;
      return;
   }

   boost::shared_ptr<AbsDecayList> prodDecayList = pbarpEnv::instance()->prodDecayList();
   std::vector<boost::shared_ptr<AbsDecay> > prodDecays = prodDecayList->getList();
   std::vector<boost::shared_ptr<AbsDecay> >::iterator it;
   bool particleFrompbarp=false;
   for(it=prodDecays.begin(); it!=prodDecays.end();++it){
      if( ( (*it)->daughter1Part()->name() == particleName) ||
	  ( (*it)->daughter2Part()->name() == particleName))
	 particleFrompbarp=true;
   }
   if(!particleFrompbarp){
      Alert << "Particles not produced from pbarp currently not supported." << endmsg;
      return;
   }

   for(Spin M1 = -J; M1 <= J; M1++){
     for(Spin M2 = -J; M2 <= J; M2++){
	Info << "Calculating Element (" << M1 << ", " << M2 << ")" << endmsg;
	calcSpinDensityMatrixElement(particleName, M1, M2);
     }
   }
}



void spinDensityHist::calcSpinDensityMatrixElement(std::string& particleName, Spin M1, Spin M2){

  std::stringstream newHistNameReal;
  newHistNameReal << particleName << "_" << M1 << "_" << M2 << "_Real";
  TH1F* newHistoReal = new TH1F(newHistNameReal.str().c_str(), newHistNameReal.str().c_str(), _nBins, -1, 1);
  std::stringstream newHistNameImag;
  newHistNameImag << particleName << "_" << M1 << "_" << M2 << "_Imag";
  TH1F* newHistoImag = new TH1F(newHistNameImag.str().c_str(), newHistNameImag.str().c_str(), _nBins, -1, 1);

  TH1F* normReal = new TH1F("normReal", "normReal", _nBins, -1, 1);
  TH1F* normImag = new TH1F("normImag", "normImag", _nBins, -1, 1);
  normReal->SetDirectory(0);
  normImag->SetDirectory(0);

  int eventsRead=0;
  std::vector<EvtData*>::iterator it;
  for(it = _dataList.begin(); it != _dataList.end(); ++it){
     complex<double> tempSpinDensity =
	boost::dynamic_pointer_cast<pbarpBaseLh>(_theLh)->calcSpinDensity(M1, M2, particleName, *it);

     fillHistogram(particleName, newHistoReal, *it, tempSpinDensity.real());
     fillHistogram(particleName, newHistoImag, *it, tempSpinDensity.imag());
     fillHistogram(particleName, normReal, *it, 1.0);
     fillHistogram(particleName, normImag, *it, 1.0);
     eventsRead++;
     if(eventsRead >= _maxEvents)
	break;
  }

  newHistoReal->Divide(normReal);
  newHistoImag->Divide(normImag);

  delete normReal;
  delete normImag;

}



void spinDensityHist::fillHistogram(std::string& particleName, TH1F* theHisto,
				    EvtData* theData, double spinDensityValue){

   boost::shared_ptr<AbsDecayList> absDecayList = pbarpEnv::instance()->absDecayList();
   std::vector<boost::shared_ptr<AbsDecay> > absDecays = absDecayList->getList();
   std::vector<boost::shared_ptr<AbsDecay> >::iterator it;

   Vector4<double> all4Vec=theData->FourVecsString["all"];
   Vector4<double> particle4Vec(0.,0.,0.,0.);

   for(it=absDecays.begin(); it!=absDecays.end(); ++it){

      if( (*it)->motherPart()->name() == particleName){
	 std::vector<Particle*> fsParticles = (*it)->finalStateParticles();
	 std::vector<Particle*>::iterator it2;
	 Vector4<double> allDaughters4Vec(0., 0., 0., 0.);

	 for(it2 = fsParticles.begin(); it2!=fsParticles.end(); ++it2){
	     allDaughters4Vec += theData->FourVecsString[(*it2)->name()];
	 }

	 particle4Vec = allDaughters4Vec;
	 particle4Vec.Boost(all4Vec);

	 break;
      }
   }

   theHisto->Fill(particle4Vec.CosTheta(), spinDensityValue * theData->evtWeight);
}
