//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
//									  //
//  This file is part of Pawian.					  //
//									  //
//  Pawian is free software: you can redistribute it and/or modify	  //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	  //
//									  //
//  Pawian is distributed in the hope that it will be useful,		  //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  //
//  GNU General Public License for more details.	      		  //
//									  //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	  //
//									  //
//************************************************************************//

// spinDensityHist class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#include <sstream>

#include "pbarpUtils/spinDensityHist.hh"
#include "pbarpUtils/pbarpEnv.hh"
#include "pbarpUtils/pbarpBaseLh.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/PwaCovMatrix.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"

#include "TH1F.h"
#include "TFile.h"


const int spinDensityHist::MAX_EVENTS = 2000;


spinDensityHist::spinDensityHist(boost::shared_ptr<AbsLh> theLh, fitParams& theFitParams) :
  _calcErrors(false)
  , _nBins(101)
  ,_theLh(theLh)
{
   _theFitParamsOriginal = &theFitParams;
   _dataList=_theLh->getMcVec();
}




spinDensityHist::~spinDensityHist(){
   _spinDensityRootFile->Write();
   _spinDensityRootFile->Close();
}



void spinDensityHist::SetCovarianceMatrix(boost::shared_ptr<PwaCovMatrix> thePwaCovMatrix){

   _calcErrors = true;
   _thePwaCovMatrix = thePwaCovMatrix;
}



void spinDensityHist::Calculate(){

   theFitParamsBaseClass.setMnUsrParams(_theMnUserParameters, *_theFitParamsOriginal, *_theFitParamsOriginal);
   
   std::stringstream spinDensityRootFileName;
   spinDensityRootFileName << "./spinDensity" << pbarpEnv::instance()->outputFileNameSuffix() << ".root";

   _spinDensityRootFile = new TFile(spinDensityRootFileName.str().c_str(), "recreate");

   std::vector<string> spinDensityParticles = pbarpEnv::instance()->spinDensityNames();

   for(auto it=spinDensityParticles.begin(); it!=spinDensityParticles.end(); ++it){
      Info << "Calculating spin density matrix for particle " << (*it) << endmsg;
      calcSpinDensityMatrix(*it);
   }
}



void spinDensityHist::calcSpinDensityMatrix(std::string& particleName){

   int J = pbarpEnv::instance()->particleTable()->particle(particleName)->J();

   if(J!=1){
      Alert << "Spin density calculation for J!=1 currently not supported." << endmsg;
      return;
   }

   boost::shared_ptr<AbsDecayList> prodDecayList = pbarpEnv::instance()->prodDecayList();
   std::vector<boost::shared_ptr<AbsDecay> > prodDecays = prodDecayList->getList();

   bool particleFrompbarp=false;
   for(auto it=prodDecays.begin(); it!=prodDecays.end();++it){
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
   
   TH1F* normHist = new TH1F("normHist", "normHist", _nBins, -1, 1);
   normHist->SetDirectory(0);
   TH1F* errReal = new TH1F("errReal", "errRead", _nBins, -1, 1);
   TH1F* errImag = new TH1F("errImag", "errImag", _nBins, -1, 1);
   errReal->SetDirectory(0);
   errImag->SetDirectory(0);
   
   int eventsRead=0;
   for(auto it = _dataList.begin(); it != _dataList.end(); ++it){

      double costheta = ParticleCosTheta(particleName, *it);

      if(newHistoReal->GetBinContent(newHistoReal->FindBin(costheta)) != 0)
	 continue;

      _theLh->updateFitParams(*_theFitParamsOriginal);
      complex<double> tempSpinDensity =
	 boost::dynamic_pointer_cast<pbarpBaseLh>(_theLh)->calcSpinDensity(M1, M2, particleName, *it);

      newHistoReal->Fill(costheta, tempSpinDensity.real() * (*it)->evtWeight);
      newHistoImag->Fill(costheta, tempSpinDensity.imag() * (*it)->evtWeight);
      normHist->Fill(costheta, (*it)->evtWeight);

      complex<double> tempSpinDensityErr(0,0);
      if(_calcErrors){
	 tempSpinDensityErr = calcSpinDensityMatrixError(particleName, M1, M2, *it, tempSpinDensity);
      }
      
      errReal->Fill(costheta, tempSpinDensityErr.real() * (*it)->evtWeight);
      errImag->Fill(costheta, tempSpinDensityErr.imag() * (*it)->evtWeight);

      eventsRead++;
      if(eventsRead >= MAX_EVENTS)
	 break;
   }

   newHistoReal->Divide(normHist);
   newHistoImag->Divide(normHist);
   errReal->Divide(normHist);
   errImag->Divide(normHist);
   
   for(int i=1; i<=_nBins; i++){
      newHistoReal->SetBinError(i, errReal->GetBinContent(i));
      newHistoImag->SetBinError(i, errImag->GetBinContent(i));
   }
   
   delete normHist;
   delete errReal;
   delete errImag;
}



complex<double> 
spinDensityHist::calcSpinDensityMatrixError(std::string& particleName, 
					    Spin M1, Spin M2, EvtData* evtData, complex<double> sdmValue){

   double stepSize = 0.0001;
   std::map< std::string, complex<double> > derivatives;

   unsigned int nPar = _theMnUserParameters.Params().size();

   for(unsigned int i=0; i<nPar; i++){

      double parOrig = _theMnUserParameters.Value(i);
      std::string parName = _theMnUserParameters.GetName(i);
      
      _theMnUserParameters.SetValue(i, parOrig + stepSize);
      
      fitParams newFitParams = *_theFitParamsOriginal;
      theFitParamsBaseClass.getFitParamVal(_theMnUserParameters.Params(), newFitParams);
      _theLh->updateFitParams(newFitParams);
      
      complex<double> tempSpinDensity  = 
	 boost::dynamic_pointer_cast<pbarpBaseLh>(_theLh)->calcSpinDensity(M1, M2, particleName, evtData);
      
      complex<double> newDerivative = (tempSpinDensity - sdmValue) / stepSize;
      derivatives[parName] = newDerivative;
      
      _theMnUserParameters.SetValue(i, parOrig);
   }

   double errorReal=0;
   double errorImag=0;

   for(unsigned int i=0; i<nPar; i++)
      for(unsigned int j=0; j<nPar; j++){
	 
	 std::string name1 = _theMnUserParameters.GetName(i);
	 std::string name2 = _theMnUserParameters.GetName(j);
	 
	 errorReal += derivatives[name1].real() * _thePwaCovMatrix->GetElement(name1, name2) * derivatives[name2].real();
	 errorImag += derivatives[name1].imag() * _thePwaCovMatrix->GetElement(name1, name2) * derivatives[name2].imag();
      }

   errorReal = sqrt(errorReal);
   errorImag = sqrt(errorImag);

   return complex<double>(errorReal, errorImag);
}



double spinDensityHist::ParticleCosTheta(std::string& particleName, EvtData* theData){

   boost::shared_ptr<AbsDecayList> absDecayList = pbarpEnv::instance()->absDecayList();
   std::vector<boost::shared_ptr<AbsDecay> > absDecays = absDecayList->getList();
   
   Vector4<double> all4Vec=theData->FourVecsString["all"];
   Vector4<double> particle4Vec(0.,0.,0.,0.);
   
   for(auto it=absDecays.begin(); it!=absDecays.end(); ++it){
      
      if( (*it)->motherPart()->name() == particleName){
	 std::vector<Particle*> fsParticles = (*it)->finalStateParticles();
	 Vector4<double> allDaughters4Vec(0., 0., 0., 0.);

	 for(auto it2 = fsParticles.begin(); it2!=fsParticles.end(); ++it2){
	    allDaughters4Vec += theData->FourVecsString[(*it2)->name()];
	 }
	 
	 particle4Vec = allDaughters4Vec;
	 particle4Vec.Boost(all4Vec);

	 break;
      }
   }
   
   return particle4Vec.CosTheta();
}
