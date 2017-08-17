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

// KMatrixDynamics class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/KMatrixDynamics.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "PwaDynamics/FVector.hh"
#include "PwaDynamics/PVectorRelBg.hh"
#include "ConfigParser/KMatrixParser.hh"
#include "ConfigParser/ParserBase.hh"
#include "PwaDynamics/KMatrixRel.hh"
#include "PwaDynamics/KMatrixRelBg.hh"
#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/KPoleBarrier.hh"
#include "PwaDynamics/PPole.hh"
#include "PwaDynamics/PPoleBarrier.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpaceFactory.hh"
#include "FitParams/AbsPawianParameters.hh"

KMatrixDynamics::KMatrixDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::string& pathToConfigParser, ChannelID channelID, std::string projectionParticleNames) :
  TMatrixDynamics(name, fsParticles, mother, pathToConfigParser, "Data", projectionParticleNames)
  ,_currentMass(1.)
  ,_pVecSuffix("")
{
  init();
  _isLdependent=true;
  ParserBase* currentParser=GlobalEnv::instance()->Channel(channelID)->parser();
  const std::vector<std::string> kMatProdSuffixes=currentParser->addKmatrixProdSuffix();
  std::map<std::string, std::string> kMatKeyProdSuffNames;

  std::vector<std::string>::const_iterator itStr;
  for ( itStr = kMatProdSuffixes.begin(); itStr != kMatProdSuffixes.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);

    std::string keyStr;
    stringStr >> keyStr;

    std::string suffStr;
    stringStr >> suffStr;
    kMatKeyProdSuffNames[keyStr]=suffStr;
  }

  std::map<std::string, std::string>::const_iterator it = kMatKeyProdSuffNames.find(_kMatName);
  if(it!=kMatKeyProdSuffNames.end()){
    _pVecSuffix=kMatKeyProdSuffNames.at(_kMatName);;
  }
}

KMatrixDynamics::~KMatrixDynamics()
{
}

complex<double> KMatrixDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){

complex<double> result(0.,0.);

  int evtNo=theData->evtNo;
  std::string currentKey="default";

  if(0!=grandmaAmp) currentKey=_massKey+grandmaAmp->absDec()->massParKey();

  if ( _cacheAmps && !_recalcMap.at(currentKey)){
     result=_cachedStringOrbMap.at(evtNo).at(currentKey).at(OrbMom);
  }
  
  else{
      theMutex.lock();
      result=_fVecMap.at(currentKey)->evalProjMatrix(theData->DoubleMassId.at(_dynId), _decProjectionIndex, OrbMom);
      //      InfoMsg << name() << "\t_decProjectionIndex: " << "\t currentKey: " << currentKey << "\t" << _decProjectionIndex << "\tresult: " << result << endmsg;
      if ( _cacheAmps){
      	_cachedStringOrbMap[evtNo][currentKey][OrbMom]=result;
      }
      theMutex.unlock();
  }
  return result;
}

void KMatrixDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  //beta factor for production
  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMap.begin(); it1!=_currentbFactorMap.end(); ++it1){
    std::string theName=it1->first;
    std::map<std::string, double>& bFactors = it1->second;
    for(unsigned int i=0; i<_poleNames.size(); ++i){ 
      std::string currentName="b"+_pVecSuffix+"_"+_poleNames.at(i);
      //     std::cout << "currentName: " << currentName << std::endl;
      std::string magName=currentName+"Mag";
      fitPar->Add(theName+magName, bFactors.at(magName) , 1.);
      //      fitPar->SetLimits(theName+magName, 0., bFactors.at(magName)+30.);

      std::string phiName=currentName+"Phi";
      fitPar->Add(theName+phiName, bFactors.at(phiName) , 0.2);
      
    }
  }
  
  //pole positions
  for(unsigned int i=0; i<_currentPoleMasses.size(); ++i){
    double valMass=_currentPoleMasses.at(i);
    double errMass=0.02;
    double minMass=valMass-5.*errMass;
    if(minMass<0.) minMass=0.;
    double maxMass=valMass+5.*errMass;
    fitPar->Add(_poleNames.at(i)+"Mass", valMass, errMass);
    fitPar->SetLimits(_poleNames.at(i)+"Mass", minMass, maxMass);
  }

    //g-factors
  for(unsigned int i=0; i<_poleNames.size(); ++i){
    std::vector<double> currentgFactorVec=_currentgFactorMap.at(i);
    for(unsigned int j=0; j<currentgFactorVec.size(); ++j){
      std::string currentName=_poleNames.at(i)+_gFactorNames.at(j)+"gFactor";
      double currentError=currentgFactorVec.at(j)/3.;
      if (currentError<1.e-5) currentError=0.01;
      fitPar->Add(currentName, currentgFactorVec.at(j), currentError);
    }
  }

    //k-matrix bg-terms
  if(_orderKMatBg>=0){
    for(unsigned int i=0; i<=fabs(_orderKMatBg); ++i){
      for(unsigned int j=0; j<_phpVecs.size(); ++j){
  	for(unsigned int k=j; k<_phpVecs.size(); ++k){
  	  std::string currentName=_bgTermNames.at(i).at(j).at(k);
  	  fitPar->Add(currentName, _currentBgTerms.at(i).at(j).at(k), fabs(_currentBgTerms.at(i).at(j).at(k))+0.3);
  	}
      }
    }
    //Adler-term
    if(_withKMatAdler){
      fitPar->Add("s0"+_kMatName, _currentAdler0, fabs(_currentAdler0)+0.2);
    }
  }

    //p-vector bgterms bg-terms
  if(_orderPVecBg>=0){
    for(unsigned int i=0; i<=fabs(_orderPVecBg); ++i){
      for(unsigned int j=0; j<_phpVecs.size(); ++j){
	std::string currentName=_bgPVecTermNames.at(i).at(j);
	fitPar->Add(currentName, _currentPVecBgTerms.at(i).at(j), fabs(_currentPVecBgTerms.at(i).at(j))+0.3);
      }
    }
  }
}

void KMatrixDynamics::fillParamNameList(){
  _paramNameList.clear();

  //beta factor for production
  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMap.begin(); it1!=_currentbFactorMap.end(); ++it1){
    std::string theName=it1->first;
    std::vector<std::string> currentNameList;
    
//    std::map<std::string, double>& bFactors = it1->second;
    for(unsigned int i=0; i<_poleNames.size(); ++i){ 
      std::string currentName="b"+_pVecSuffix+"_"+_poleNames.at(i);
      //     std::cout << "currentName: " << currentName << std::endl;
      std::string magName=currentName+"Mag";
      _paramNameList.push_back(theName+magName);
      currentNameList.push_back(theName+magName);
  
      std::string phiName=currentName+"Phi";
      _paramNameList.push_back(theName+phiName);
      currentNameList.push_back(theName+phiName);
    }
    _paramNameListMap[theName]=currentNameList;
  }

  std::map<std::string, std::vector<std::string> >::iterator itNameMap;  
  //pole positions
  for(unsigned int i=0; i<_currentPoleMasses.size(); ++i){
    _paramNameList.push_back(_poleNames.at(i)+"Mass");
    for(itNameMap=_paramNameListMap.begin(); itNameMap!=_paramNameListMap.end(); ++itNameMap){
      itNameMap->second.push_back(_poleNames.at(i)+"Mass");
    }
  }

    //g-factors
  for(unsigned int i=0; i<_poleNames.size(); ++i){
    std::vector<double> currentgFactorVec=_currentgFactorMap.at(i);
    for(unsigned int j=0; j<currentgFactorVec.size(); ++j){
      std::string currentName=_poleNames.at(i)+_gFactorNames.at(j)+"gFactor";
      _paramNameList.push_back(currentName);
      for(itNameMap=_paramNameListMap.begin(); itNameMap!=_paramNameListMap.end(); ++itNameMap){
	itNameMap->second.push_back(currentName);
      }
    }
  }

    //k-matrix bg-terms
  if(_orderKMatBg>=0){
    for(unsigned int i=0; i<=fabs(_orderKMatBg); ++i){
      for(unsigned int j=0; j<_phpVecs.size(); ++j){
  	for(unsigned int k=j; k<_phpVecs.size(); ++k){
  	  std::string currentName=_bgTermNames.at(i).at(j).at(k);
	  _paramNameList.push_back(currentName);
          for(itNameMap=_paramNameListMap.begin(); itNameMap!=_paramNameListMap.end(); ++itNameMap){
	    itNameMap->second.push_back(currentName);
	  }
  	}
      }
    }
    //Adler-term
    if(_withKMatAdler){
      _paramNameList.push_back("s0"+_kMatName);
      for(itNameMap=_paramNameListMap.begin(); itNameMap!=_paramNameListMap.end(); ++itNameMap){
	itNameMap->second.push_back("s0"+_kMatName);
      }
    }
  }

 //p-vector bg-terms
  if(_orderPVecBg>=0){
    for(unsigned int i=0; i<=fabs(_orderPVecBg); ++i){
      for(unsigned int j=0; j<_phpVecs.size(); ++j){
  	  std::string currentName=_bgPVecTermNames.at(i).at(j);
	  _paramNameList.push_back(currentName);
          for(itNameMap=_paramNameListMap.begin(); itNameMap!=_paramNameListMap.end(); ++itNameMap){
	    itNameMap->second.push_back(currentName);
  	}
      }
    }
  }
}

bool KMatrixDynamics::checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld){

  std::map<std::string, std::vector<std::string> >::iterator itMap;
  for(itMap= _paramNameListMap.begin(); itMap !=_paramNameListMap.end(); ++itMap){
    _recalcMap.at(itMap->first)=false;
    std::vector<std::string>::iterator itStr;
    for (itStr=itMap->second.begin(); itStr!=itMap->second.end(); ++itStr){
      std::string currentParamName=*itStr;
      if(!CheckDoubleEquality(fitParNew->Value(currentParamName), fitParOld->Value(currentParamName))){
        _recalcMap.at(itMap->first)=true;
        continue;
      }
    }
  }

  return AbsParamHandler::checkRecalculation(fitParNew, fitParOld);
}

void KMatrixDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  //beta factor for production  
  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMap.begin(); it1!=_currentbFactorMap.end(); ++it1){
    std::map<std::string, double>::iterator it2;
    std::map<std::string, double>& bFactors=it1->second;
    for(it2=bFactors.begin(); it2!=bFactors.end(); ++it2){
      it2->second=fitPar->Value(it1->first+it2->first);
    }

    std::shared_ptr<PVectorRel> currentPVec=_pVecMap.at(it1->first);

    for(unsigned int i=0; i<_poleNames.size(); ++i){ 
      std::string currentName="b"+_pVecSuffix+"_"+_poleNames.at(i);
      complex<double> currentbFactor=fabs(bFactors.at(currentName+"Mag"))*complex<double>(cos(bFactors.at(currentName+"Phi")), sin(bFactors.at(currentName+"Phi")));
      currentPVec->updateBeta(i, currentbFactor);
    }
  }

  //pole positions
  for(unsigned int i=0; i<_currentPoleMasses.size(); ++i){
    std::string currentPoleName=_poleNames.at(i)+"Mass";
    double currentPoleMass=fitPar->Value(currentPoleName);
    _currentPoleMasses.at(i)=currentPoleMass;
    _kPoles.at(i)->updatePoleMass(currentPoleMass);

    std::map<std::string, std::shared_ptr<PVectorRel> >::iterator itPVec;
    for(itPVec=_pVecMap.begin(); itPVec!=_pVecMap.end(); ++itPVec){
      itPVec->second->updatePoleMass(i, currentPoleMass);
    }
  }

  //g-factors
  for(unsigned int i=0; i<_poleNames.size(); ++i){
    std::vector<double>& currentgFactorVec=_currentgFactorMap.at(i);
    for(unsigned int j=0; j<currentgFactorVec.size(); ++j){
      std::string currentName=_poleNames.at(i)+_gFactorNames.at(j)+"gFactor";
      currentgFactorVec.at(j)=fitPar->Value(currentName);
    }
    _kPoles.at(i)->updategFactors(currentgFactorVec);
    
    std::map<std::string, std::shared_ptr<PVectorRel> >::iterator itPVec;
    for(itPVec=_pVecMap.begin(); itPVec!=_pVecMap.end(); ++itPVec){
      itPVec->second->updategFactors(i, currentgFactorVec);
    }
  }

  //k-matrix bg-terms
  if(_orderKMatBg>=0){
    for(unsigned int i=0; i<=fabs(_orderKMatBg); ++i){
      for(unsigned int j=0; j<_phpVecs.size(); ++j){
	for(unsigned int k=j; k<_phpVecs.size(); ++k){
	  double newVal=fitPar->Value(_bgTermNames.at(i).at(j).at(k));
	  _currentBgTerms.at(i).at(j).at(k)=newVal;
	  _kMatr->updateBgTerms(i,j,k,newVal);
	}
      }
    }

    //Adler-term
    if(_withKMatAdler){
      _currentAdler0=fitPar->Value("s0"+_kMatName);
      _kMatr->updates0Adler(_currentAdler0);
    }
  }

  //p-vector bg-terms
  if(_orderPVecBg>=0){
    for(unsigned int i=0; i<=fabs(_orderPVecBg); ++i){
      for(unsigned int j=0; j<_phpVecs.size(); ++j){
	  double newVal=fitPar->Value(_bgPVecTermNames.at(i).at(j));
	  _currentPVecBgTerms.at(i).at(j)=newVal;
	  std::map<std::string, std::shared_ptr<PVectorRel> >::iterator itPVec;
	  for(itPVec=_pVecMap.begin(); itPVec!=_pVecMap.end(); ++itPVec){
	    itPVec->second->updateBgTerms(i, j, newVal);
	  }
      }
    }
  }
}


void KMatrixDynamics::addGrandMa(std::shared_ptr<AbsDecay> theDec){
  if(0==theDec){
     Alert << "Can not add AbsXdecAmp; 0 pointer!!!" << endmsg;
     exit(1);
  }
  
  std::string theName=_massKey+theDec->massParKey();

  std::cout << "addGrandMa:\t" << theName << std::endl;

  std::map<std::string, std::shared_ptr<FVector> >::iterator it = _fVecMap.find(theName);
  
  if (it != _fVecMap.end()) return;

  std::shared_ptr<PVectorRel> currentPVector=makeNewPVec();
  _pVecMap[theName]=currentPVector;

  std::vector< std::string>::iterator poleNameIt;
  for (poleNameIt=_poleNames.begin(); poleNameIt!=_poleNames.end(); ++poleNameIt){  
    std::string currentName="b"+_pVecSuffix+"_"+(*poleNameIt);
    _currentbFactorMap[theName][currentName+"Mag"]=1.;
    _currentbFactorMap[theName][currentName+"Phi"]=0.;
  } 

  std::map<std::string, double>& bFactors = _currentbFactorMap[theName];
  for(unsigned int i=0; i<_poleNames.size(); ++i){ 
    std::string currentName="b"+_pVecSuffix+"_"+_poleNames[i];
    std::cout << "currentName: " << currentName << std::endl;
    complex<double> currentbFactor=bFactors[currentName+"Mag"]*complex<double>(cos(bFactors[currentName+"Phi"]), sin(bFactors[currentName+"Phi"]));
    currentPVector->updateBeta(i, currentbFactor);
  }

  if(_orderPVecBg>=0) {
    _currentPVecBgTerms.resize(_orderPVecBg+1);
    _bgPVecTermNames.resize(_orderPVecBg+1);
    for(unsigned int i=0; i<= fabs(_orderPVecBg); ++i){
      _currentPVecBgTerms.at(i).resize(_phpVecs.size());
      _bgPVecTermNames.at(i).resize(_phpVecs.size());
      for(unsigned int j=0; j<_phpVecs.size(); ++j){
	_currentPVecBgTerms.at(i).at(j)=0.;

	  std::stringstream keyOrderStrStr;
	  keyOrderStrStr << i << j;
	  std::string keyOrder=keyOrderStrStr.str();
	  std::string currentName="bgPVec"+_pVecSuffix+keyOrder+_kMatName;
	  _bgPVecTermNames.at(i).at(j)=currentName;
      }
    }
  }


  std::shared_ptr<FVector> currentFVector=std::shared_ptr<FVector>(new FVector(_kMatr, currentPVector));
  _fVecMap[theName]=currentFVector;
  _recalcMap[theName]=true;
}

const std::string& KMatrixDynamics::grandMaKey(AbsXdecAmp* grandmaAmp){
  if(0==grandmaAmp) return _grandmaKey;
  return grandmaAmp->absDec()->massParKey();
}

const unsigned short KMatrixDynamics::grandMaId(AbsXdecAmp* grandmaAmp){
  if(0==grandmaAmp) return _grandmaId;
  return grandmaAmp->absDec()->massParId(); 
}


void KMatrixDynamics::init(){
  _orderPVecBg=_kMatrixParser->orderPVecBg();
}


std::shared_ptr<PVectorRel> KMatrixDynamics::makeNewPVec(){

  vector<std::shared_ptr<PPole> > thePpoles;
  complex<double> defaultBeta(1.,0.); 
     
   vector<std::shared_ptr<KPole> >::iterator it;
   for (it=_kPoles.begin(); it!=_kPoles.end(); ++it){ 
     std::vector<double> currentGFactors=(*it)->gFactors();
     std::shared_ptr<PPole> currentPPole;
     if (_kMatrixParser->useBarrierFactors()) currentPPole=std::shared_ptr<PPole>(new PPoleBarrier(defaultBeta, currentGFactors, (*it)->poleMass(), _phpVecs, _kMatrixParser->orbitalMom(), _kMatrixParser->useTruncatedBarrierFactors()));
     else currentPPole=std::shared_ptr<PPole>(new PPole(defaultBeta, currentGFactors, (*it)->poleMass())); 
     thePpoles.push_back(currentPPole);     
   }

   std::shared_ptr<PVectorRel> thePVector;
   if(_orderPVecBg<0) thePVector=std::shared_ptr<PVectorRel>(new PVectorRel(thePpoles, _phpVecs));
   else thePVector=std::shared_ptr<PVectorRel>(new PVectorRelBg(thePpoles, _phpVecs, _orderPVecBg)); 
   return thePVector;
}
