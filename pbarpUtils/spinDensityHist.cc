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
#include "pbarpUtils/pbarpBaseLh.hh"
#include "pbarpUtils/PbarpChannelEnv.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "FitParams/PwaCovMatrix.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "FitParams/ParamFactory.hh"
#include "Utils/IdStringMapRegistry.hh"

#include "TH1F.h"
#include "TFile.h"


const int spinDensityHist::MAX_EVENTS = 10000;


spinDensityHist::spinDensityHist(std::shared_ptr<AbsLh> theLh, std::shared_ptr<AbsPawianParameters> theFitParams) :
   _calcErrors(false)
  , _nBins(101)
  ,_theLh(theLh)
   ,_theFitParamsOriginal(theFitParams)
   ,_theParameters(ParamFactory::instance()->getParametersPointer("Pawian"))
{
  //   _theFitParamsOriginal = &theFitParams;
   _dataList=_theLh->getMcVec();
}




spinDensityHist::~spinDensityHist(){
   _spinDensityRootFile->Write();
   _spinDensityRootFile->Close();
}



void spinDensityHist::SetCovarianceMatrix(std::shared_ptr<PwaCovMatrix> thePwaCovMatrix){

   _calcErrors = true;
   _thePwaCovMatrix = thePwaCovMatrix;
}



void spinDensityHist::Calculate(){
  _theParameters=std::shared_ptr<AbsPawianParameters>(_theFitParamsOriginal->Clone());
  //  theFitParColBaseClass.setAbsPawianParams(_theParameters, *_theFitParamsOriginal, *_theFitParamsOriginal);
   
   std::stringstream spinDensityRootFileName;
   spinDensityRootFileName << "./spinDensity" << GlobalEnv::instance()->outputFileNameSuffix() << ".root";

   _spinDensityRootFile = new TFile(spinDensityRootFileName.str().c_str(), "recreate");

   std::vector<string> spinDensityParticles = 
      static_pointer_cast<PbarpChannelEnv>(GlobalEnv::instance()->PbarpChannel())->spinDensityNames();

   for(auto it=spinDensityParticles.begin(); it!=spinDensityParticles.end(); ++it){
      Info << "Calculating spin density matrix for particle " << (*it) << endmsg;
      calcSpinDensityMatrix(*it);
   }
}



void spinDensityHist::calcSpinDensityMatrix(std::string& particleName){

   int J = GlobalEnv::instance()->particleTable()->particle(particleName)->J();

   if(J<=0){
       Alert << "Cannot calculate spin density matrix for spin 0 particles" << endmsg;
       return;
   }

   if(J>1){
       Warning << "Spin density calculation for J!=1 experimental." << endmsg;
   }

   std::shared_ptr<AbsDecayList> prodDecayList = 
      static_pointer_cast<PbarpChannelEnv>(GlobalEnv::instance()->PbarpChannel())->prodDecayList();
   std::vector<std::shared_ptr<AbsDecay> > prodDecays = prodDecayList->getList();

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

   calcSpinDensityMatrixNorm(particleName, J);
   for(Spin M1 = -J; M1 <= J; ++M1){
      for(Spin M2 = -J; M2 <= J; ++M2){
	Info << "Calculating Element (" << M1 << ", " << M2 << ")" << endmsg;
	calcSpinDensityMatrixElement(particleName, M1, M2, J);
      }
   }
}



void spinDensityHist::calcSpinDensityMatrixNorm(std::string& particleName, int J){

   std::stringstream newHistNameNorm;
   newHistNameNorm << particleName << "_norm";
   TH1F* normHist = new TH1F(newHistNameNorm.str().c_str(), newHistNameNorm.str().c_str(), _nBins, -1, 1);

   TH1F* normnormHist = new TH1F("normnormHist", "normnormHist", _nBins, -1, 1);
   normnormHist->SetDirectory(0);

   std::vector<short> entries (_nBins);
   std::fill(entries.begin(), entries.end(), 0);

   for(auto it = _dataList.begin(); it != _dataList.end(); ++it){

      double costheta = ParticleCosTheta(particleName, *it);
      short histIndex = normHist->FindBin(costheta);

      if(entries.at(histIndex-1) >= 2)
	 continue;

      entries.at(histIndex-1)++;

      _theLh->updateFitParams(_theFitParamsOriginal);
      double tempSpinDensityNorm =
	 std::dynamic_pointer_cast<pbarpBaseLh>(_theLh)->calcSpinDensityNorm(particleName, *it, J);

      normHist->Fill(costheta, tempSpinDensityNorm * (*it)->evtWeight);
      normnormHist->Fill(costheta, (*it)->evtWeight);
   }

   normHist->Divide(normnormHist);
   delete normnormHist;
}



void spinDensityHist::calcSpinDensityMatrixElement(std::string& particleName, Spin M1, Spin M2, int J){

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

   std::vector<short> entries (_nBins);
   std::fill(entries.begin(), entries.end(), 0);

   for(auto it = _dataList.begin(); it != _dataList.end(); ++it){

      double costheta = ParticleCosTheta(particleName, *it);
      short histIndex = newHistoReal->FindBin(costheta);

      if(entries.at(histIndex-1) >= 2)
	 continue;

      entries.at(histIndex-1)++;

      _theLh->updateFitParams(_theFitParamsOriginal);
      complex<double> tempSpinDensity =
	 std::dynamic_pointer_cast<pbarpBaseLh>(_theLh)->calcSpinDensity(M1, M2, particleName, *it, J);

      newHistoReal->Fill(costheta, tempSpinDensity.real() * (*it)->evtWeight);
      newHistoImag->Fill(costheta, tempSpinDensity.imag() * (*it)->evtWeight);
      normHist->Fill(costheta, (*it)->evtWeight);

      complex<double> tempSpinDensityErr(0,0);
      if(_calcErrors){
	 tempSpinDensityErr = calcSpinDensityMatrixError(particleName, M1, M2, *it, tempSpinDensity, J);
      }

      errReal->Fill(costheta, tempSpinDensityErr.real() * (*it)->evtWeight);
      errImag->Fill(costheta, tempSpinDensityErr.imag() * (*it)->evtWeight);
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
					    Spin M1, Spin M2, EvtData* evtData, complex<double> sdmValue, int J){

   std::map< std::string, complex<double> > derivatives;

   unsigned int nPar = _theParameters->Params().size();

   for(unsigned int i=0; i<nPar; i++){

      double parOrig = _theParameters->Value(i);
      std::string parName = _theParameters->GetName(i);
      double stepSize = sqrt(_thePwaCovMatrix->GetElement(parName, parName)) * 0.01;

      if(AbsParamHandler::CheckDoubleEquality(stepSize, 0)){
	 derivatives[parName] = std::complex<double>(0,0);
	 continue;
      }

      _theParameters->SetValue(i, parOrig + stepSize);
      
      //      std::shared_ptr<AbsPawianParameters> newFitParams = std::shared_ptr<AbsPawianParameters>(_theFitParamsOriginal->Clone());
      //      theFitParColBaseClass.getFitParamVal(_theParameters->Params(), newFitParams);
      //      _theLh->updateFitParams(newFitParams);
      
      _theLh->updateFitParams(_theParameters);
      
      complex<double> tempSpinDensity  = 
	 std::dynamic_pointer_cast<pbarpBaseLh>(_theLh)->calcSpinDensity(M1, M2, particleName, evtData, J);
      
      complex<double> newDerivative = (tempSpinDensity - sdmValue) / stepSize;
      derivatives[parName] = newDerivative;
      
      _theParameters->SetValue(i, parOrig);
   }

   double errorReal=0;
   double errorImag=0;

   for(unsigned int i=0; i<nPar; i++)
      for(unsigned int j=0; j<nPar; j++){
	 
	 std::string name1 = _theParameters->GetName(i);
	 std::string name2 = _theParameters->GetName(j);
	 
	 errorReal += derivatives[name1].real() * _thePwaCovMatrix->GetElement(name1, name2) * derivatives[name2].real();
	 errorImag += derivatives[name1].imag() * _thePwaCovMatrix->GetElement(name1, name2) * derivatives[name2].imag();
      }

   errorReal = sqrt(errorReal);
   errorImag = sqrt(errorImag);

   return complex<double>(errorReal, errorImag);
}



double spinDensityHist::ParticleCosTheta(std::string& particleName, EvtData* theData){

   std::shared_ptr<AbsDecayList> absDecayList = static_pointer_cast<PbarpChannelEnv>(GlobalEnv::instance()->PbarpChannel())->absDecayList();
   std::vector<std::shared_ptr<AbsDecay> > absDecays = absDecayList->getList();
   
   std::string stringAll="all";
   Vector4<double> all4Vec=theData->FourVecsId.at(IdStringMapRegistry::instance()->stringId(stringAll));

   Vector4<double> particle4Vec(0.,0.,0.,0.);
   
   for(auto it=absDecays.begin(); it!=absDecays.end(); ++it){
      
      if( (*it)->motherPart()->name() == particleName){
	 std::vector<Particle*> fsParticles = (*it)->finalStateParticles();
	 Vector4<double> allDaughters4Vec(0., 0., 0., 0.);

	 for(auto it2 = fsParticles.begin(); it2!=fsParticles.end(); ++it2){
	    std::string currentName=(*it2)->name();
	    allDaughters4Vec += theData->FourVecsId.at(IdStringMapRegistry::instance()->stringId(currentName));
	 }
	 
	 particle4Vec = allDaughters4Vec;
	 particle4Vec.Boost(all4Vec);

	 break;
      }
   }
   
   return particle4Vec.CosTheta();
}


