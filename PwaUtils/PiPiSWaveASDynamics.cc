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
#include <thread>

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
#include "FitParams/AbsPawianParameters.hh"

PiPiSWaveASDynamics::PiPiSWaveASDynamics(std::string& name, std::vector<Particle*>& fsParticles, 
					 Particle* mother, ParticleTable* thePdtTable) :
  AbsDynamics(name, fsParticles, mother)
  ,_pdtTable(thePdtTable)
  ,_projectionIndex(projectionIndex(fsParticles))
  , _KmatrixPiPiS(new KMatrixPiPiS())
  , _phpVec(_KmatrixPiPiS->phaseSpaceVec())
  ,_kPoles(_KmatrixPiPiS->kpoles())
{
  DebugMsg << "projection index is " << _projectionIndex << endmsg;
  _isLdependent=false;
}

PiPiSWaveASDynamics::~PiPiSWaveASDynamics()
{
}

complex<double> PiPiSWaveASDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){

complex<double> result(0.,0.);

  int evtNo=theData->evtNo;
  std::string currentKey="default";

  if(0!=grandmaAmp) currentKey=_massKey+grandmaAmp->absDec()->massParKey();

  if ( _cacheAmps){
    if (!_recalcMap.at(currentKey)) result= _cachedStringMap.at(evtNo).at(currentKey); 
      else{ 
	theMutex.lock();
	result=_fVecMap.at(currentKey)->evalProjMatrix(theData->DoubleMassId.at(_dynId), _projectionIndex);
	_cachedStringMap[evtNo][currentKey]=result; 
	theMutex.unlock();
      }
  }
  else if (!_cacheAmps) result=_fVecMap.at(currentKey)->evalProjMatrix(theData->DoubleMassId.at(_dynId), _projectionIndex);
  
  return result;
}

void PiPiSWaveASDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  //beta factor for production
  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMagMap.begin(); it1!=_currentbFactorMagMap.end(); ++it1){
    std::string theName=it1->first;
    std::map<std::string, double>& bMagFactors = it1->second;

    std::map<std::string, double>::iterator it2;    
    for(it2=bMagFactors.begin(); it2!=bMagFactors.end(); ++it2){
      std::string magName=theName+it2->first;
      fitPar->Add(magName, it2->second , 1.);
    }

    std::map<std::string, double>& bPhiFactors = _currentbFactorPhiMap.at(it1->first);
    for(it2=bPhiFactors.begin(); it2!=bPhiFactors.end(); ++it2){
      std::string phiName=theName+it2->first;
      fitPar->Add(phiName, 0. , 0.2);
    }

    //fProd factors
    std::map<std::string, double>& fMagProds=_currentfProdMagMap.at(it1->first);
    for(it2=fMagProds.begin(); it2!=fMagProds.end(); ++it2){
      std::string magName=theName+it2->first;
      fitPar->Add(magName, it2->second , 1.);
    }

    std::map<std::string, double>& fPhiProds=_currentfProdPhiMap.at(it1->first);
    for(it2=fPhiProds.begin(); it2!=fPhiProds.end(); ++it2){
      std::string phiName=theName+it2->first;
      fitPar->Add(phiName, 0., 0.2);
    }

    //s0 factors
    std::string s0Name=theName+"S0";
    fitPar->Add(s0Name, _currentS0Map.at(it1->first) , 0.2);
  }
}

void PiPiSWaveASDynamics::fillParamNameList(){
  //beta factor for production
  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMagMap.begin(); it1!=_currentbFactorMagMap.end(); ++it1){
    std::string theName=it1->first;

    std::vector<std::string> currentNameList;
    
    std::map<std::string, double>& bMagFactors = it1->second;

    std::map<std::string, double>::iterator it2;    
    for(it2=bMagFactors.begin(); it2!=bMagFactors.end(); ++it2){
      std::string magName=theName+it2->first;
      _paramNameList.push_back(magName);
      currentNameList.push_back(magName);
    }

    std::map<std::string, double>& bPhiFactors = _currentbFactorPhiMap.at(it1->first);
    for(it2=bPhiFactors.begin(); it2!=bPhiFactors.end(); ++it2){
      std::string phiName=theName+it2->first;
      _paramNameList.push_back(phiName);
      currentNameList.push_back(phiName);
    }

    //fProd factors
    std::map<std::string, double>& fMagProds=_currentfProdMagMap.at(it1->first);
    for(it2=fMagProds.begin(); it2!=fMagProds.end(); ++it2){
      std::string magName=theName+it2->first;
      _paramNameList.push_back(magName);
      currentNameList.push_back(magName);
    }

    std::map<std::string, double>& fPhiProds=_currentfProdPhiMap.at(it1->first);
    for(it2=fPhiProds.begin(); it2!=fPhiProds.end(); ++it2){
      std::string phiName=theName+it2->first;
      _paramNameList.push_back(phiName);
      currentNameList.push_back(phiName);
    }

    //s0 factors
    std::string s0Name=theName+"S0";
    _paramNameList.push_back(s0Name);
    currentNameList.push_back(s0Name);
    _paramNameListMap[theName]=currentNameList;  
  }

  InfoMsg << "***** list of parameter names for PiPiSWaveASDynamics ******" << endmsg;
  std::vector<std::string>::iterator itstr;
  for (itstr=_paramNameList.begin(); itstr!=_paramNameList.end(); ++itstr){
    InfoMsg << (*itstr) << endmsg;
  } 
}

bool PiPiSWaveASDynamics::checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, 
					     std::shared_ptr<AbsPawianParameters> fitParOld){

  std::map<std::string, std::vector<std::string> >::iterator itMap;
  for(itMap= _paramNameListMap.begin(); itMap !=_paramNameListMap.end(); ++itMap){
    _recalcMap[itMap->first]=false;
    std::vector<std::string>::iterator itStr;
    for (itStr=itMap->second.begin(); itStr!=itMap->second.end(); ++itStr){
      std::string currentParamName=*itStr;
      if(!CheckDoubleEquality(fitParNew->Value(currentParamName), fitParOld->Value(currentParamName))){
	_recalcMap[itMap->first]=true;
	continue;
      }
    }
  }

 return AbsParamHandler::checkRecalculation(fitParNew, fitParOld);
}

void PiPiSWaveASDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  //beta factor for production
  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMagMap.begin(); it1!=_currentbFactorMagMap.end(); ++it1){
    std::string theName=it1->first;
    std::map<std::string, double>& bMagFactors = it1->second;

    std::map<std::string, double>::iterator it2;    
    for(it2=bMagFactors.begin(); it2!=bMagFactors.end(); ++it2){
      std::string magName=theName+it2->first;
      it2->second = fabs(fitPar->Value(magName));
    }

    std::map<std::string, double>& bPhiFactors = _currentbFactorPhiMap.at(it1->first);
    for(it2=bPhiFactors.begin(); it2!=bPhiFactors.end(); ++it2){
      std::string phiName=theName+it2->first;
      it2->second = fitPar->Value(phiName);
    }

    //fProd factors
    std::map<std::string, double>& fMagProds=_currentfProdMagMap.at(it1->first);
    for(it2=fMagProds.begin(); it2!=fMagProds.end(); ++it2){
      std::string magName=theName+it2->first;
      it2->second = fabs(fitPar->Value(magName));
    }

    std::map<std::string, double>& fPhiProds=_currentfProdPhiMap.at(it1->first);
    for(it2=fPhiProds.begin(); it2!=fPhiProds.end(); ++it2){
      std::string phiName=theName+it2->first;
      it2->second = fitPar->Value(phiName);
    }

    //s0 factors
    std::string s0Name=theName+"S0";
    _currentS0Map[it1->first]=fitPar->Value(s0Name);

    std::shared_ptr<PVectorSlowCorRel> currentPVec=_pVecMap[it1->first];

    //update _pipiSFVec
    complex<double> b_pole1=bMagFactors.at("b_pole1Mag")*complex<double>(cos(bPhiFactors.at("b_pole1Phi")), 
									 sin(bPhiFactors.at("b_pole1Phi")));
    complex<double> b_pole2=bMagFactors.at("b_pole2Mag")*complex<double>(cos(bPhiFactors.at("b_pole2Phi")), 
									 sin(bPhiFactors.at("b_pole2Phi"))); 
    complex<double> b_pole3=bMagFactors.at("b_pole3Mag")*complex<double>(cos(bPhiFactors.at("b_pole3Phi")), 
									 sin(bPhiFactors.at("b_pole3Phi")));
    complex<double> b_pole4=bMagFactors.at("b_pole4Mag")*complex<double>(cos(bPhiFactors.at("b_pole4Phi")), 
									 sin(bPhiFactors.at("b_pole4Phi")));
    complex<double> b_pole5=bMagFactors.at("b_pole5Mag")*complex<double>(cos(bPhiFactors.at("b_pole5Phi")), 
									 sin(bPhiFactors.at("b_pole5Phi")));
    currentPVec->updateBeta(0, b_pole1);
    currentPVec->updateBeta(1, b_pole2);
    currentPVec->updateBeta(2, b_pole3);
    currentPVec->updateBeta(3, b_pole4);
    currentPVec->updateBeta(4, b_pole5);

    complex<double> fProdPiPi =
      fMagProds.at("fprod_PiPiMag")*complex<double>(cos(fPhiProds.at("fprod_PiPiPhi")), 
						    sin(fPhiProds.at("fprod_PiPiPhi")));
    complex<double> fProdKK = 
      fMagProds.at("fprod_KKMag")*complex<double>(cos(fPhiProds.at("fprod_KKPhi")), 
						  sin(fPhiProds.at("fprod_KKPhi")));
    complex<double> fProd4Pi =
      fMagProds.at("fprod_4PiMag")*complex<double>(cos(fPhiProds.at("fprod_4PiPhi")), 
						   sin(fPhiProds.at("fprod_4PiPhi")));
    complex<double> fProdEtaEta =
      fMagProds.at("fprod_EtaEtaMag")*complex<double>(cos(fPhiProds.at("fprod_EtaEtaPhi")), 
						      sin(fPhiProds.at("fprod_EtaEtaPhi")));
    complex<double> fProdEtaEtap =
      fMagProds.at("fprod_EtaEtapMag")*complex<double>(cos(fPhiProds.at("fprod_EtaEtapPhi")), 
						       sin(fPhiProds.at("fprod_EtaEtapPhi")));
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

  _currentbFactorMagMap[theName]["b_pole1Mag"]=1.;
  _currentbFactorPhiMap[theName]["b_pole1Phi"]=0.;
  _currentbFactorMagMap[theName]["b_pole2Mag"]=1.;
  _currentbFactorPhiMap[theName]["b_pole2Phi"]=0.;
  _currentbFactorMagMap[theName]["b_pole3Mag"]=1.;
  _currentbFactorPhiMap[theName]["b_pole3Phi"]=0.;
  _currentbFactorMagMap[theName]["b_pole4Mag"]=1.;
  _currentbFactorPhiMap[theName]["b_pole4Phi"]=0.;
  _currentbFactorMagMap[theName]["b_pole5Mag"]=1.;
  _currentbFactorPhiMap[theName]["b_pole5Phi"]=0.;

  std::map<std::string, double>& bMagFactors = _currentbFactorMagMap[theName];
  std::map<std::string, double>& bPhiFactors = _currentbFactorPhiMap[theName];
  complex<double> b_pole1=bMagFactors["b_pole1Mag"]*complex<double>(cos(bPhiFactors["b_pole1Phi"]), 
								    sin(bPhiFactors["b_pole1Phi"]));
  complex<double> b_pole2=bMagFactors["b_pole2Mag"]*complex<double>(cos(bPhiFactors["b_pole2Phi"]), 
								    sin(bPhiFactors["b_pole2Phi"])); 
  complex<double> b_pole3=bMagFactors["b_pole3Mag"]*complex<double>(cos(bPhiFactors["b_pole3Phi"]), 
								    sin(bPhiFactors["b_pole3Phi"]));
  complex<double> b_pole4=bMagFactors["b_pole4Mag"]*complex<double>(cos(bPhiFactors["b_pole4Phi"]), 
								    sin(bPhiFactors["b_pole4Phi"]));
  complex<double> b_pole5=bMagFactors["b_pole5Mag"]*complex<double>(cos(bPhiFactors["b_pole5Phi"]), 
								    sin(bPhiFactors["b_pole5Phi"]));
  currentPVector->updateBeta(0, b_pole1);
  currentPVector->updateBeta(1, b_pole2);
  currentPVector->updateBeta(2, b_pole3);
  currentPVector->updateBeta(3, b_pole4);
  currentPVector->updateBeta(4, b_pole5);

  _currentfProdMagMap[theName]["fprod_PiPiMag"]=1.;
  _currentfProdPhiMap[theName]["fprod_PiPiPhi"]=0.;
  _currentfProdMagMap[theName]["fprod_KKMag"]=1.;
  _currentfProdPhiMap[theName]["fprod_KKPhi"]=0.;
  _currentfProdMagMap[theName]["fprod_4PiMag"]=1.;
  _currentfProdPhiMap[theName]["fprod_4PiPhi"]=0.;
  _currentfProdMagMap[theName]["fprod_EtaEtaMag"]=1.;
  _currentfProdPhiMap[theName]["fprod_EtaEtaPhi"]=0.;
  _currentfProdMagMap[theName]["fprod_EtaEtapMag"]=1.;
  _currentfProdPhiMap[theName]["fprod_EtaEtapPhi"]=0.;

  std::map<std::string, double>& fMagProds=_currentfProdMagMap[theName];
  std::map<std::string, double>& fPhiProds=_currentfProdPhiMap[theName];
  complex<double> fProdPiPi=fMagProds["fprod_PiPiMag"]*complex<double>(cos(fPhiProds["fprod_PiPiPhi"]), 
								       sin(fPhiProds["fprod_PiPiPhi"]));
  complex<double> fProdKK=fMagProds["fprod_KKMag"]*complex<double>(cos(fPhiProds["fprod_KKPhi"]), 
								   sin(fPhiProds["fprod_KKPhi"]));
  complex<double> fProd4Pi=fMagProds["fprod_4PiMag"]*complex<double>(cos(fPhiProds["fprod_4PiPhi"]), 
								     sin(fPhiProds["fprod_4PiPhi"]));
  complex<double> fProdEtaEta=fMagProds["fprod_EtaEtaMag"]*complex<double>(cos(fPhiProds["fprod_EtaEtaPhi"]), 
									   sin(fPhiProds["fprod_EtaEtaPhi"]));
  complex<double> fProdEtaEtap = 
    fMagProds["fprod_EtaEtapMag"]*complex<double>(cos(fPhiProds["fprod_EtaEtapPhi"]), 
						  sin(fPhiProds["fprod_EtaEtapPhi"]));

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

const unsigned short PiPiSWaveASDynamics::grandMaId(AbsXdecAmp* grandmaAmp){
  if(0==grandmaAmp) return _grandmaId;
  return grandmaAmp->absDec()->massParId(); 
}

int PiPiSWaveASDynamics::projectionIndex(std::vector<Particle*>& fsParticles){

  int result=-1;

  if ( fsParticles.size() != 2 && fsParticles.size() != 4){
    Alert << "PiPiSWave: decay to " << FunctionUtils::particleListName(fsParticles) 
	  << " is not supported!!!" << endmsg;
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
    Alert << "PiPiSWave: decay to " << FunctionUtils::particleListName(fsParticles) 
	  << " is not supported!!!" <<endmsg;
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
