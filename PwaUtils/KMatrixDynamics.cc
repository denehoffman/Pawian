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
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "PwaDynamics/FVector.hh"
#include "ConfigParser/KMatrixParser.hh"
#include "PwaDynamics/KMatrixRel.hh"
#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/KPoleBarrier.hh"
#include "PwaDynamics/PPole.hh"
#include "PwaDynamics/PPoleBarrier.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpaceIsobar.hh"

KMatrixDynamics::KMatrixDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::string& pathToConfigParser) :
  AbsDynamics(name, fsParticles, mother)
  ,_kMatrixParser(new KMatrixParser(pathToConfigParser))
{
  init();
}

KMatrixDynamics::~KMatrixDynamics()
{
}

complex<double> KMatrixDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){

complex<double> result(0.,0.);

  int evtNo=theData->evtNo;
  std::string currentKey="default";

  if(0!=grandmaAmp) currentKey=_massKey+grandmaAmp->absDec()->massParKey();

  if ( _cacheAmps && !_recalcMap[currentKey]){
    result=_cachedStringMap[evtNo][currentKey];
  }
  else{
      theMutex.lock();
      result=_fVecMap[currentKey]->evalProjMatrix(theData->FourVecsString[_dynKey].M(), _projectionIndex);
      if ( _cacheAmps){
  	_cachedStringMap[evtNo][currentKey]=result;
      }
      theMutex.unlock();
  }
  
  return result;
}

void  KMatrixDynamics::getDefaultParams(fitParams& fitVal, fitParams& fitErr){

  //beta factor for production  
  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMap.begin(); it1!=_currentbFactorMap.end(); ++it1){
    
    std::map<std::string, double>::iterator it2;
    
    std::map<std::string, double>& bFactors=it1->second;
    for(it2=bFactors.begin(); it2!=bFactors.end(); ++it2){
      fitVal.otherParams[it1->first+it2->first]=it2->second;
      fitErr.otherParams[it1->first+it2->first]=1.;
    }
  }

  //pole positions
  std::vector<double >::iterator itPoleVec;
  for(unsigned int i=0; i<_currentPoleMasses.size(); ++i){
    //  for(itPoleVec=_currentPoleMasses.begin(); itPoleVec!=_currentPoleMasses.end(); ++itPoleVec){
    fitVal.Masses[_poleNames.at(i)]=_currentPoleMasses.at(i);
    fitErr.Masses[_poleNames.at(i)]=0.02;
  }

  //g-factors
  for(unsigned int i=0; i<_poleNames.size(); ++i){
    std::vector<double> currentgFactorVec=_currentgFactorMap.at(i);
    for(unsigned int j=0; j<currentgFactorVec.size(); ++j){
      std::string currentName=_poleNames.at(i)+_gFactorNames.at(j)+"gFactor";
      fitVal.gFactors[currentName]=currentgFactorVec.at(j);
      fitErr.gFactors[currentName]=currentgFactorVec.at(j);
    }
  }
 

  }

bool KMatrixDynamics::checkRecalculation(fitParams& theParamVal){
  _recalculate=false;

  //beta factor for production  
  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMap.begin(); it1!=_currentbFactorMap.end(); ++it1){
    _recalcMap[it1->first]=false;    
    std::map<std::string, double>::iterator it2;
    
    std::map<std::string, double>& bFactors=it1->second;
    for(it2=bFactors.begin(); it2!=bFactors.end(); ++it2){
      if (!CheckDoubleEquality(it2->second, theParamVal.otherParams.at(it1->first+it2->first))){
	_recalcMap[it1->first]=true;
      }
    }
  }

  //pole positions
  std::vector<double >::iterator itPoleVec;
  for(unsigned int i=0; i<_currentPoleMasses.size(); ++i){
    if (!CheckDoubleEquality( _currentPoleMasses.at(i), theParamVal.Masses[_poleNames.at(i)])){
      _recalculate=true;
    }
  }

  //g-factors
  for(unsigned int i=0; i<_poleNames.size(); ++i){
    std::vector<double> currentgFactorVec=_currentgFactorMap.at(i);
    for(unsigned int j=0; j<currentgFactorVec.size(); ++j){
      std::string currentName=_poleNames.at(i)+_gFactorNames.at(j)+"gFactor";
      if (!CheckDoubleEquality( currentgFactorVec.at(j), theParamVal.gFactors.at(currentName))){
	_recalculate=true;
      }
    }
  }

  if (_recalculate){
    std::map<std::string, bool >::iterator itRecalcMap;
    for(itRecalcMap=_recalcMap.begin(); itRecalcMap!=_recalcMap.end(); ++itRecalcMap){
      itRecalcMap->second=true;
    } 
  }

  return _recalculate;
}

void KMatrixDynamics::updateFitParams(fitParams& theParamVal){

  //beta factor for production  
  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMap.begin(); it1!=_currentbFactorMap.end(); ++it1){
    std::map<std::string, double>::iterator it2;
    std::map<std::string, double>& bFactors=it1->second;
    for(it2=bFactors.begin(); it2!=bFactors.end(); ++it2){
      it2->second=theParamVal.otherParams.at(it1->first+it2->first);
    }

    std::shared_ptr<PVectorRel> currentPVec=_pVecMap[it1->first];

    for(unsigned int i=0; i<_poleNames.size(); ++i){ 
      std::string currentName="b_"+_poleNames.at(i);
      //      std::cout << "currentName: " << currentName << std::endl;
      complex<double> currentbFactor=bFactors.at(currentName+"Mag")*complex<double>(cos(bFactors.at(currentName+"Phi")), sin(bFactors.at(currentName+"Phi")));
      currentPVec->updateBeta(i, currentbFactor);
    }
  }

  //pole positions
  std::vector<double >::iterator itPoleVec;
  for(unsigned int i=0; i<_currentPoleMasses.size(); ++i){
    double currentPoleMass=theParamVal.Masses[_poleNames.at(i)];
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
      currentgFactorVec.at(j)=theParamVal.gFactors.at(currentName);
    }
    _kPoles.at(i)->updategFactors(currentgFactorVec);
    
    std::map<std::string, std::shared_ptr<PVectorRel> >::iterator itPVec;
    for(itPVec=_pVecMap.begin(); itPVec!=_pVecMap.end(); ++itPVec){
      itPVec->second->updategFactors(i, currentgFactorVec);
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
    std::string currentName="b_"+(*poleNameIt);
    _currentbFactorMap[theName][currentName+"Mag"]=1.;
    _currentbFactorMap[theName][currentName+"Phi"]=0.;
  } 

  std::map<std::string, double>& bFactors = _currentbFactorMap[theName];
  for(unsigned int i=0; i<_poleNames.size(); ++i){ 
    std::string currentName="b_"+_poleNames[i];
    std::cout << "currentName: " << currentName << std::endl;
    complex<double> currentbFactor=bFactors[currentName+"Mag"]*complex<double>(cos(bFactors[currentName+"Phi"]), sin(bFactors[currentName+"Phi"]));
    currentPVector->updateBeta(i, currentbFactor);
  }

  std::shared_ptr<FVector> currentFVector=std::shared_ptr<FVector>(new FVector(_kMatr, currentPVector));
  _fVecMap[theName]=currentFVector;
  _recalcMap[theName]=true;
}

const std::string& KMatrixDynamics::grandMaKey(AbsXdecAmp* grandmaAmp){
  if(0==grandmaAmp) return _grandmaKey;
  //  return grandmaAmp->absDec()->dynKey();
  return grandmaAmp->absDec()->massParKey();
}



void KMatrixDynamics::init(){
  
  std::vector<std::string> poleNameAndMassVecs=_kMatrixParser->poles();
  std::vector<std::string>::iterator itString;
  for (itString=poleNameAndMassVecs.begin(); itString!=poleNameAndMassVecs.end(); ++itString){
    std::istringstream poleIString(*itString);
    std::string currentPoleName;
    std::string currentPoleMassStr;
    poleIString >> currentPoleName >> currentPoleMassStr;

    std::istringstream currentPoleMassiStr(currentPoleMassStr);
    double currentValue;
    if(!(currentPoleMassiStr >> currentValue)){
      Alert << "cannot convert " << currentPoleMassStr << " to a double value" << endmsg;
      exit(0);
    }

    _currentPoleMasses.push_back(currentValue);
    _poleNames.push_back(currentPoleName);
  }

  if(_currentPoleMasses.size()!= _kMatrixParser->noOfPoles()){
    Alert << "number of poles != number of pole masses" << endmsg;
    exit(0);
  }

  const std::vector<std::string> gFacStringVec=_kMatrixParser->gFactors();
  std::vector<std::string>::const_iterator itStrConst;
  for(itStrConst=gFacStringVec.begin(); itStrConst!=gFacStringVec.end(); ++itStrConst){
    std::istringstream particles(*itStrConst);
    std::string firstParticleName;
    std::string secondParticleName;
    particles >> firstParticleName >> secondParticleName;

    Particle* firstParticle = GlobalEnv::instance()->particleTable()->particle(firstParticleName);
    Particle* secondParticle = GlobalEnv::instance()->particleTable()->particle(secondParticleName);
    if(0==firstParticle || 0==secondParticle){
      Alert << "particle with name: " << firstParticleName <<" or " << secondParticleName << " doesn't exist in pdg-table" << endmsg;
      exit(0);
    }
    std::shared_ptr<AbsPhaseSpace> currentPhp(new PhaseSpaceIsobar(firstParticle->mass(), secondParticle->mass()));
    _phpVecs.push_back(currentPhp);

    std::string gFactorKey=firstParticleName+secondParticleName;    
    _gFactorNames.push_back(gFactorKey);   

    for (int i=0; i<int(_kMatrixParser->noOfPoles()); ++i){
      std::string currentPoleName=_poleNames[i];
      std::string currentgValueStr;
      if(!(particles >> currentgValueStr)){
	Alert << "g-factors for pole " << currentPoleName << " does not exist!" << endmsg;
	exit(0);
      }
      std::istringstream currentgValueiStr(currentgValueStr);
      double currentGValue;
      if (!(currentgValueiStr >> currentGValue)){
	Alert << "cannot convert " << currentgValueStr << " to a double value" << endmsg;
	exit(0);
      }
      std::string gFactorKey=firstParticleName+secondParticleName;
      _currentgFactorMap[i].push_back(currentGValue);
    }
  }


  std::map<int, std::vector<double> >::iterator itgFac;
  for (itgFac=_currentgFactorMap.begin(); itgFac!=_currentgFactorMap.end(); ++itgFac){
    std::vector<double> currentgVector=itgFac->second;
    std::shared_ptr<KPole> currentPole;
    if (_kMatrixParser->useBarrierFactors()) currentPole=std::shared_ptr<KPole>(new KPoleBarrier(currentgVector, _currentPoleMasses.at(itgFac->first), _phpVecs, _kMatrixParser->orbitalMom()));
    else currentPole=std::shared_ptr<KPole>(new KPole(currentgVector, _currentPoleMasses.at(itgFac->first)));
    _kPoles.push_back(currentPole);
  }

  _kMatr=std::shared_ptr<KMatrixRel>(new KMatrixRel(_kPoles,_phpVecs ));

  const std::string porjectionParticleNames=_kMatrixParser->projection();
  std::istringstream projParticles(porjectionParticleNames);
  std::string firstProjParticleName;
  std::string secondProjParticleName;
  projParticles >> firstProjParticleName >> secondProjParticleName;
  std::string projKey=firstProjParticleName+secondProjParticleName;    

  bool found=false;
  for(unsigned int i=0; i<_gFactorNames.size();++i){
    if(projKey==_gFactorNames[i]){
      _projectionIndex=i;
      found=true;
    }
  }
  if (!found){
    Alert << "projection index for key " << projKey << " not found" << endmsg;
    exit(0);
  }  
}


std::shared_ptr<PVectorRel> KMatrixDynamics::makeNewPVec(){

  vector<std::shared_ptr<PPole> > thePpoles;
  complex<double> defaultBeta(1.,0.); 
     
   vector<std::shared_ptr<KPole> >::iterator it;
   for (it=_kPoles.begin(); it!=_kPoles.end(); ++it){ 
     std::vector<double> currentGFactors=(*it)->gFactors();
     std::shared_ptr<PPole> currentPPole;
     if (_kMatrixParser->useBarrierFactors()) currentPPole=std::shared_ptr<PPole>(new PPoleBarrier(defaultBeta, currentGFactors, (*it)->poleMass(), _phpVecs, _kMatrixParser->orbitalMom()));
     else currentPPole=std::shared_ptr<PPole>(new PPole(defaultBeta, currentGFactors, (*it)->poleMass())); 
     thePpoles.push_back(currentPPole);     
   }

   std::shared_ptr<PVectorRel> thePVector(new PVectorRel(thePpoles, _phpVecs));
   return thePVector;
}
