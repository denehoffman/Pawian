//************************************************************************//
//									  //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// OmnesDynamics class definition file. -*- C++ -*-
// Copyright 2020 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/OmnesDynamics.hh"
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
#include "PwaDynamics/KMatrixOmnes.hh"
#include "PwaDynamics/KPoleOmnes.hh"
#include "PwaDynamics/PPoleOmnes.hh"
#include "PwaDynamics/FVectorOmnes.hh"
#include "PwaDynamics/PVectorBgOmnes.hh"
#include "PwaDynamics/PhaseSpaceList.hh"
#include "qft++Extension/AbsOmnesMatrix.hh"
#include "qft++Extension/OmnesMatrixFactory.hh"
#include "FitParams/AbsPawianParameters.hh"

OmnesDynamics::OmnesDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::string& pathToConfigParser, ChannelID channelID, std::string projectionParticleNames) :
  KMatrixDynamics(name, fsParticles, mother, pathToConfigParser, channelID, projectionParticleNames)
{
  init();
  if(_orderPVecBg<0){
    Alert << "OmnesDynamics requires P-vector background terms; _orderPVecBg = "
	  << _orderPVecBg << " is not allowed" << endmsg;
    exit(1);
  }
}

OmnesDynamics::~OmnesDynamics()
{
}

complex<double> OmnesDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){

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

void OmnesDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
    //beta factor for production
  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMap.begin(); it1!=_currentbFactorMap.end(); ++it1){
    std::string theName=it1->first;
    std::map<std::string, double>& bFactors = it1->second;
    for(unsigned int i=0; i<_poleNames.size(); ++i){ 
      std::string currentName="b"+_pVecSuffix+"_"+_poleNames.at(i);
      std::string magName=currentName+"Mag";
      fitPar->Add(theName+magName, bFactors.at(magName) , 1.);
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

void OmnesDynamics::fillParamNameList(){
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

void OmnesDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
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

std::string OmnesDynamics::addOneGrandMa(std::string theName){
  std::shared_ptr<PVectorBgOmnes> currentPVector=makeNewPVecOmnes();
  _pVecMap[theName]=currentPVector;
  std::cout << "theName: " << theName << std::endl;

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
	  if(_useParticleNameForPVecBg) currentName="bgPVec"+_pVecSuffix+keyOrder+theName;
	  _bgPVecTermNames.at(i).at(j)=currentName;
      }
    }
  }

  //std::shared_ptr<PhaseSpaceList> phpList(new PhaseSpaceList(_phpVecs));
  std::string omnesMatrType=_kMatrixParser->omnesMatrixType();
  std::shared_ptr<AbsOmnesMatrix> onmesMatrPointer=OmnesMatrixFactory::instance()->getOmnesMatrixPointer(omnesMatrType);
  std::shared_ptr<PhaseSpaceList> phpList(new PhaseSpaceList(_phpVecs, onmesMatrPointer->omnesMatrix(), onmesMatrPointer->selfEnergyMatrix()));
  std::shared_ptr<FVector> currentFVector=std::shared_ptr<FVector>(new FVectorOmnes(_kMatrOmnes, currentPVector, phpList));
  
  _fVecMap[theName]=currentFVector;
  _recalcMap[theName]=true;

  return theName;
}


std::shared_ptr<PVectorBgOmnes> OmnesDynamics::makeNewPVecOmnes(){

  vector<std::shared_ptr<PPole> > thePpoles;
  complex<double> defaultBeta(1.,0.); 
     
   vector<std::shared_ptr<KPoleOmnes> >::iterator it;
   for (it=_kPolesOmnes.begin(); it!=_kPolesOmnes.end(); ++it){ 
     std::vector<double> currentGFactors=(*it)->gFactors();
     std::shared_ptr<PPoleOmnes> currentPPole(new PPoleOmnes(defaultBeta, currentGFactors, (*it)->poleMass())); 
     thePpoles.push_back(currentPPole);     
   }

   std::shared_ptr<PVectorBgOmnes> thePVector(new PVectorBgOmnes(thePpoles, _phpVecs, _orderPVecBg)); 
   return thePVector;
}

void OmnesDynamics::init(){
  std::map<int, std::vector<double> >::iterator itgFac;
  for (itgFac=_currentgFactorMap.begin(); itgFac!=_currentgFactorMap.end(); ++itgFac){
    std::vector<double> currentgVector=itgFac->second;
    std::shared_ptr<KPoleOmnes> currentPole;
    currentPole=std::shared_ptr<KPoleOmnes>(new KPoleOmnes(currentgVector, _currentPoleMasses.at(itgFac->first)));
    _kPolesOmnes.push_back(currentPole);
  }

  _kMatrOmnes=std::shared_ptr<KMatrixOmnes>(new KMatrixOmnes(_kPoles,_phpVecs ));

}

