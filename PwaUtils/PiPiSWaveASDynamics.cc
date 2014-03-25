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

// PiPiSWaveASDynamics class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/PiPiSWaveASDynamics.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Utils/FunctionUtils.hh"
#include "PwaDynamics/FVector.hh"
#include "PwaDynamics/KMatrixPiPiS.hh" 
#include "PwaDynamics/PVectorSlowCorRel.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/PPole.hh"
#include "Particle/ParticleTable.hh"

PiPiSWaveASDynamics::PiPiSWaveASDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, ParticleTable* thePdtTable) :
  AbsDynamics(name, fsParticles, mother)
  ,_pdtTable(thePdtTable)
  ,_projectionIndex(projectionIndex(fsParticles))
  , _KmatrixPiPiS(new KMatrixPiPiS())
  , _phpVec(_KmatrixPiPiS->phaseSpaceVec())
  ,_kPoles(_KmatrixPiPiS->kpoles())
{
  DebugMsg << "projection index is " << _projectionIndex << endmsg;
}

PiPiSWaveASDynamics::~PiPiSWaveASDynamics()
{
}

complex<double> PiPiSWaveASDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){

complex<double> result(0.,0.);

  int evtNo=theData->evtNo;
  std::string currentKey="default";

  if(0!=grandmaAmp) currentKey=_massKey+grandmaAmp->absDec()->massParKey();
 
  bool readFromCachedMap=false;

  if( _cacheAmps){
    if(_recalcMap[currentKey]){
      bool currentEvtAlreadyCached=false;
      
      std::map<int, std::map<std::string, bool > >::iterator itAlreadyCached=_alreadyCached.find(evtNo);
      if( itAlreadyCached != _alreadyCached.end()){
	std::map<std::string, bool >::iterator itAlreadyCached2= itAlreadyCached->second.find(currentKey);
	if( itAlreadyCached2 != itAlreadyCached->second.end()){
	  currentEvtAlreadyCached=itAlreadyCached2->second;     
	}
	else{
	  theMutex.lock();
	  _alreadyCached[evtNo][currentKey]=false;
	  theMutex.unlock();
	}   
      }
      else{ 
	theMutex.lock();
	_alreadyCached[evtNo][currentKey]=false;
	theMutex.unlock();
      } 
  
      if(currentEvtAlreadyCached) readFromCachedMap=true;
    }
    else readFromCachedMap=true; 
  }

  if ( readFromCachedMap){
    result=_cachedStringMap.at(evtNo).at(currentKey);
  }  
  else{
      theMutex.lock();
      result=_fVecMap[currentKey]->evalProjMatrix(theData->FourVecsString[_dynKey].M(), _projectionIndex);
      if ( _cacheAmps){
	_cachedStringMap[evtNo][currentKey]=result;
	_alreadyCached.at(evtNo).at(currentKey)=true;
      }
      theMutex.unlock();
  }
  //  Info << "result: " << result << endmsg;  
  return result;
}

void  PiPiSWaveASDynamics::getDefaultParams(fitParams& fitVal, fitParams& fitErr){
  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMap.begin(); it1!=_currentbFactorMap.end(); ++it1){
    
    std::map<std::string, double>::iterator it2;
    
    std::map<std::string, double>& bFactors=it1->second;
    for(it2=bFactors.begin(); it2!=bFactors.end(); ++it2){
      fitVal.otherParams[it1->first+it2->first]=it2->second;
      fitErr.otherParams[it1->first+it2->first]=1.;
    }

    std::map<std::string, double>& fProds=_currentfProdMap[it1->first];
    for(it2=fProds.begin(); it2!=fProds.end(); ++it2){
      fitVal.otherParams[it1->first+it2->first]=it2->second;
      fitErr.otherParams[it1->first+it2->first]=1.; 
    }

    fitVal.otherParams[it1->first+"S0_PosNeg"]=_currentS0Map[it1->first];
    fitErr.otherParams[it1->first+"S0_PosNeg"]=4.;
  }
}

bool PiPiSWaveASDynamics::checkRecalculation(fitParams& theParamVal){
  _recalculate=false;
  std::map<int, std::map<std::string, bool > >::iterator itAlreadyCached;
  for(itAlreadyCached=_alreadyCached.begin(); itAlreadyCached!=_alreadyCached.end(); ++itAlreadyCached){
    std::map<std::string, bool >::iterator itAlreadyCached2;
    for(itAlreadyCached2=itAlreadyCached->second.begin(); itAlreadyCached2!=itAlreadyCached->second.end(); ++itAlreadyCached2){
      itAlreadyCached2->second=false;     
    }
  }  


  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMap.begin(); it1!=_currentbFactorMap.end(); ++it1){

    _recalcMap[it1->first]=false;    
    
     std::map<std::string, double>::iterator it2;
    
    std::map<std::string, double>& bFactors=it1->second;
    for(it2=bFactors.begin(); it2!=bFactors.end(); ++it2){
       if (!CheckDoubleEquality(it2->second, theParamVal.otherParams[it1->first+it2->first])){
	_recalculate=true;
	_recalcMap[it1->first]=true;
	continue;
      }
    }

   std::map<std::string, double>& fProds=_currentfProdMap[it1->first];
   for(it2=fProds.begin(); it2!=fProds.end(); ++it2){
     if(_recalcMap[it1->first]==true) continue;
     if (!CheckDoubleEquality(it2->second, theParamVal.otherParams[it1->first+it2->first])){
       _recalculate=true;
       _recalcMap[it1->first]=true;
     }
   }

   if (!CheckDoubleEquality(_currentS0Map[it1->first], theParamVal.otherParams[it1->first+"S0_PosNeg"])){
     _recalculate=true;
     _recalcMap[it1->first]=true;
   }
   
  }
  
  return _recalculate;
}

void PiPiSWaveASDynamics::updateFitParams(fitParams& theParamVal){

  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMap.begin(); it1!=_currentbFactorMap.end(); ++it1){

    std::map<std::string, double>::iterator it2;

    std::map<std::string, double>& bFactors = it1->second;
    for(it2=bFactors.begin(); it2!=bFactors.end(); ++it2){
      it2->second = theParamVal.otherParams[it1->first+it2->first];
    }

    std::map<std::string, double>& fProds=_currentfProdMap[it1->first];
    for(it2=fProds.begin(); it2!=fProds.end(); ++it2){
      it2->second = theParamVal.otherParams[it1->first+it2->first];
    }

    _currentS0Map[it1->first]=theParamVal.otherParams[it1->first+"S0_PosNeg"];

    // std::shared_ptr<FVectorPiPiS> currentFVec=_fVecMap[it1->first];
    std::shared_ptr<PVectorSlowCorRel> currentPVec=_pVecMap[it1->first];

    //update _pipiSFVec
    complex<double> b_pole1=bFactors["b_pole1Mag"]*complex<double>(cos(bFactors["b_pole1Phi"]), sin(bFactors["b_pole1Phi"]));
    complex<double> b_pole2=bFactors["b_pole2Mag"]*complex<double>(cos(bFactors["b_pole2Phi"]), sin(bFactors["b_pole2Phi"])); 
    complex<double> b_pole3=bFactors["b_pole3Mag"]*complex<double>(cos(bFactors["b_pole3Phi"]), sin(bFactors["b_pole3Phi"]));
    complex<double> b_pole4=bFactors["b_pole4Mag"]*complex<double>(cos(bFactors["b_pole4Phi"]), sin(bFactors["b_pole4Phi"]));
    complex<double> b_pole5=bFactors["b_pole5Mag"]*complex<double>(cos(bFactors["b_pole5Phi"]), sin(bFactors["b_pole5Phi"]));
    currentPVec->updateBeta(0, b_pole1);
    currentPVec->updateBeta(1, b_pole2);
    currentPVec->updateBeta(2, b_pole3);
    currentPVec->updateBeta(3, b_pole4);
    currentPVec->updateBeta(4, b_pole5);

    complex<double> fProdPiPi=fProds["fprod_PiPiMag"]*complex<double>(cos(fProds["fprod_PiPiPhi"]), sin(fProds["fprod_PiPiPhi"]));
    complex<double> fProdKK=fProds["fprod_KKMag"]*complex<double>(cos(fProds["fprod_KKPhi"]), sin(fProds["fprod_KKPhi"]));
    complex<double> fProd4Pi=fProds["fprod_4PiMag"]*complex<double>(cos(fProds["fprod_4PiPhi"]), sin(fProds["fprod_4PiPhi"]));
    complex<double> fProdEtaEta=fProds["fprod_EtaEtaMag"]*complex<double>(cos(fProds["fprod_EtaEtaPhi"]), sin(fProds["fprod_EtaEtaPhi"]));
    complex<double> fProdEtaEtap=fProds["fprod_EtaEtapMag"]*complex<double>(cos(fProds["fprod_EtaEtapPhi"]), sin(fProds["fprod_EtaEtapPhi"]));
    currentPVec->updateFprod(0, fProdPiPi);
    currentPVec->updateFprod(1, fProdKK);
    currentPVec->updateFprod(2, fProd4Pi);
    currentPVec->updateFprod(3, fProdEtaEta);
    currentPVec->updateFprod(4, fProdEtaEtap);
    
    currentPVec->updateS0prod(_currentS0Map[it1->first]);
   }
}

void PiPiSWaveASDynamics::addGrandMa(std::shared_ptr<AbsDecay> theDec){
  if(0==theDec){
    Alert << "Can not add AbsXdecAmp; 0 pointer!!!" << endmsg;
    exit(1);
  }
  
  std::string theName=_massKey+theDec->massParKey();

  std::cout << "addGrandMa:\t" << theName << std::endl;

  std::map<std::string, std::shared_ptr<FVector> >::iterator it = _fVecMap.find(theName);
  
  if (it != _fVecMap.end()) return;

  std::shared_ptr<PVectorSlowCorRel> currentPVector=makeNewPVec();
  std::shared_ptr<FVector> currentFVector(new FVector(_KmatrixPiPiS, currentPVector));

  _currentbFactorMap[theName]["b_pole1Mag"]=1.;
  _currentbFactorMap[theName]["b_pole1Phi"]=0.;
  _currentbFactorMap[theName]["b_pole2Mag"]=1.;
  _currentbFactorMap[theName]["b_pole2Phi"]=0.;
  _currentbFactorMap[theName]["b_pole3Mag"]=1.;
  _currentbFactorMap[theName]["b_pole3Phi"]=0.;
  _currentbFactorMap[theName]["b_pole4Mag"]=1.;
  _currentbFactorMap[theName]["b_pole4Phi"]=0.;
  _currentbFactorMap[theName]["b_pole5Mag"]=1.;
  _currentbFactorMap[theName]["b_pole5Phi"]=0.;

  std::map<std::string, double>& bFactors = _currentbFactorMap[theName];
  complex<double> b_pole1=bFactors["b_pole1Mag"]*complex<double>(cos(bFactors["b_pole1Phi"]), sin(bFactors["b_pole1Phi"]));
  complex<double> b_pole2=bFactors["b_pole2Mag"]*complex<double>(cos(bFactors["b_pole2Phi"]), sin(bFactors["b_pole2Phi"])); 
  complex<double> b_pole3=bFactors["b_pole3Mag"]*complex<double>(cos(bFactors["b_pole3Phi"]), sin(bFactors["b_pole3Phi"]));
  complex<double> b_pole4=bFactors["b_pole4Mag"]*complex<double>(cos(bFactors["b_pole4Phi"]), sin(bFactors["b_pole4Phi"]));
  complex<double> b_pole5=bFactors["b_pole5Mag"]*complex<double>(cos(bFactors["b_pole5Phi"]), sin(bFactors["b_pole5Phi"]));
  currentPVector->updateBeta(0, b_pole1);
  currentPVector->updateBeta(1, b_pole2);
  currentPVector->updateBeta(2, b_pole3);
  currentPVector->updateBeta(3, b_pole4);
  currentPVector->updateBeta(4, b_pole5);

  _currentfProdMap[theName]["fprod_PiPiMag"]=1.;
  _currentfProdMap[theName]["fprod_PiPiPhi"]=0.;
  _currentfProdMap[theName]["fprod_KKMag"]=1.;
  _currentfProdMap[theName]["fprod_KKPhi"]=0.;
  _currentfProdMap[theName]["fprod_4PiMag"]=1.;
  _currentfProdMap[theName]["fprod_4PiPhi"]=0.;
  _currentfProdMap[theName]["fprod_EtaEtaMag"]=1.;
  _currentfProdMap[theName]["fprod_EtaEtaPhi"]=0.;
  _currentfProdMap[theName]["fprod_EtaEtapMag"]=1.;
  _currentfProdMap[theName]["fprod_EtaEtapPhi"]=0.;

  std::map<std::string, double>& fProds=_currentfProdMap[theName];
  complex<double> fProdPiPi=fProds["fprod_PiPiMag"]*complex<double>(cos(fProds["fprod_PiPiPhi"]), sin(fProds["fprod_PiPiPhi"]));
  complex<double> fProdKK=fProds["fprod_KKMag"]*complex<double>(cos(fProds["fprod_KKPhi"]), sin(fProds["fprod_KKPhi"]));
  complex<double> fProd4Pi=fProds["fprod_4PiMag"]*complex<double>(cos(fProds["fprod_4PiPhi"]), sin(fProds["fprod_4PiPhi"]));
  complex<double> fProdEtaEta=fProds["fprod_EtaEtaMag"]*complex<double>(cos(fProds["fprod_EtaEtaPhi"]), sin(fProds["fprod_EtaEtaPhi"]));
  complex<double> fProdEtaEtap=fProds["fprod_EtaEtapMag"]*complex<double>(cos(fProds["fprod_EtaEtapPhi"]), sin(fProds["fprod_EtaEtapPhi"]));

  currentPVector->updateFprod(0, fProdPiPi);
  currentPVector->updateFprod(1, fProdKK);
  currentPVector->updateFprod(2, fProd4Pi);
  currentPVector->updateFprod(3, fProdEtaEta);
  currentPVector->updateFprod(4, fProdEtaEtap);

  _currentS0Map[theName]=-3.;
  currentPVector->updateS0prod(_currentS0Map[theName]);

  _pVecMap[theName]=currentPVector;
  _fVecMap[theName]=currentFVector;
  _recalcMap[theName]=true;
}

const std::string& PiPiSWaveASDynamics::grandMaKey(AbsXdecAmp* grandmaAmp){
  if(0==grandmaAmp) return _grandmaKey;
  return grandmaAmp->absDec()->massParKey();
}

int PiPiSWaveASDynamics::projectionIndex(std::vector<Particle*>& fsParticles){

  int result=-1;

  if ( fsParticles.size() != 2 && fsParticles.size() != 4){
    Alert << "PiPiSWave: decay to " << FunctionUtils::particleListName(fsParticles) << " is not supported!!!" << endmsg;
    exit(0); 
  }
  if ( fsParticles.size()==4){
    result=2; //4 pi
  }
  // here  fsParticles.size()==2
  
  else{
    Particle* Piplus=_pdtTable->particle("pion+");
    Particle* Piminus=_pdtTable->particle("pion-");
    Particle* Pi0=_pdtTable->particle("pion0");
    Particle* Kplus=_pdtTable->particle("K+");
    Particle* Kminus=_pdtTable->particle("K-");
    Particle* Eta=_pdtTable->particle("eta");
    Particle* Etaprime=_pdtTable->particle("eta'");

    //first: check pi pi channel    
    if( *(fsParticles[0])==*Piplus && *(fsParticles[1])==*Piminus) result=0;
    else if(*(fsParticles[1])==*Piplus && *(fsParticles[0])==*Piminus) result=0;
    else if(*(fsParticles[0])==*Pi0 && *(fsParticles[1])==*Pi0) result=0;
    
    //second: check KK channel
    else if(*(fsParticles[0])==*Kplus && *(fsParticles[1])==*Kminus) result=1;
    else if (*(fsParticles[1])==*Kplus && *(fsParticles[0])==*Kminus) result=1;
    
    //3: check eta eta channel
    else if(*(fsParticles[0])==*Eta && *(fsParticles[1])==*Eta) result=3;
    
    //3: check eta etaprime channel
    else if(*(fsParticles[0])==*Eta && *(fsParticles[1])==*Etaprime) result=4;
    else if(*(fsParticles[1])==*Eta && *(fsParticles[0])==*Etaprime) result=4;
  }
  
  if(result==-1){
    Alert << "PiPiSWave: decay to " << FunctionUtils::particleListName(fsParticles) << " is not supported!!!" <<endmsg;
    exit(0);
  }
  
  return result;
}

std::shared_ptr<PVectorSlowCorRel> PiPiSWaveASDynamics::makeNewPVec(){

  vector<std::shared_ptr<PPole> > thePpoles;
  complex<double> defaultBeta(1.,0.); 
     
   vector<std::shared_ptr<KPole> >::iterator it;
   for (it=_kPoles.begin(); it!=_kPoles.end(); ++it){ 
     std::vector<double> currentGFactors=(*it)->gFactors();
     std::shared_ptr<PPole> currentPPole(new PPole(defaultBeta, currentGFactors, (*it)->poleMass()));
     thePpoles.push_back(currentPPole);     
   }

   std::vector<complex <double> > fProdVec;
   
   for (int i=0; i<int(_phpVec.size()); ++i){
     complex<double> currentVal(1.0,0.);
     fProdVec.push_back(currentVal);
   }
   
   double s0Prod=-0.0737;
   std::shared_ptr<PVectorSlowCorRel> thePVector(new PVectorSlowCorRel(thePpoles, _phpVec, fProdVec, s0Prod));
   return thePVector;
}
