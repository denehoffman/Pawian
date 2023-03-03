//************************************************************************//
//									  //
//  Copyright 2023 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// FixedKMatrixDynamics class definition file. -*- C++ -*-
// Copyright 2023 Bertram Kopf

#include "PwaUtils/FixedKMatrixDynamics.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "Utils/PawianIOUtils.hh"
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
#include "FitParams/AbsPawianParamStreamer.hh"

FixedKMatrixDynamics::FixedKMatrixDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::string& pathToConfigParser, ChannelID channelID, std::string projectionParticleNames) :
  KMatrixDynamics(name, fsParticles, mother, pathToConfigParser, channelID, projectionParticleNames)
{
  std::string pathToFixedKMatParams=PawianIOUtils::getFileName(GlobalEnv::instance()->KMatrixStorePath(), _kMatrixParser->fixedParamFile()); 
  AbsPawianParamStreamer thePawianStreamer(pathToFixedKMatParams);
  _kMatFixParams = thePawianStreamer.paramList();
  InfoMsg << "the fix K-matrix parameters are: " << endmsg;
  _kMatFixParams->print(std::cout);
  TMatrixDynamics::updateFitParams(_kMatFixParams);

  for(unsigned int i=0; i<_currentPoleMasses.size(); ++i){
      std::map<std::string, std::shared_ptr<PVectorRel> >::iterator itPVec;
    for(itPVec=_pVecMap.begin(); itPVec!=_pVecMap.end(); ++itPVec){
      itPVec->second->updatePoleMass(i, _currentPoleMasses.at(i));
    }
  }

}

FixedKMatrixDynamics::~FixedKMatrixDynamics()
{
}

void FixedKMatrixDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
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

void FixedKMatrixDynamics::fillParamNameList(){
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

 //p-vector bg-terms
  std::map<std::string, std::vector<std::string> >::iterator itNameMap;
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

void FixedKMatrixDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
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

