//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
//          	   - Ruhr-Universit??t Bochum 				  //
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

// TMatrixDynamics class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/TMatrixDynamics.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/EvtDataScatteringList.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "ConfigParser/KMatrixParser.hh"
#include "PwaDynamics/KMatrixRel.hh"
#include "PwaDynamics/KMatrixRelBg.hh"
#include "PwaDynamics/TMatrixRel.hh"
#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/KPoleBarrier.hh"
#include "PwaDynamics/PPole.hh"
#include "PwaDynamics/PPoleBarrier.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpaceFactory.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "Utils/IdStringMapRegistry.hh"
#include "Utils/PawianConstants.hh"

TMatrixDynamics::TMatrixDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::string& pathToConfigParser, std::string dataType) :
  AbsDynamics(name, fsParticles, mother)
  ,_kMatName("")
  ,_prodProjectionIndex(0)
  ,_decProjectionIndex(0)
  , _orderKMatBg(-1)
  ,_withKMatAdler(false)
  ,_dataTypeID(0)
  ,_prodIsNotDecChannel(false)
  ,_currentAdler0(0.)
  ,_kMatrixParser(new KMatrixParser(pathToConfigParser))
{
  if(dataType=="Elasticity") _dataTypeID=1;
  else if(dataType=="Phase") _dataTypeID=2;
  else if(dataType=="ArgandUnits") _dataTypeID=3;
  else{
    Alert << "production formalism/data type with the name" << dataType << " is not supported for pi pi scattering fits! \n It is working for: Elasticity, ArgandUnits or Phase!!!" << endmsg;
    exit(1); 
  }

  init();
  if(_prodProjectionIndex != _decProjectionIndex) _prodIsNotDecChannel=true;
  _isLdependent=true;
}

TMatrixDynamics::~TMatrixDynamics()
{
}

complex<double> TMatrixDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){

  vector<std::shared_ptr<AbsPhaseSpace> > thePhpVecs=_tMatr->kMatrix()->phaseSpaceVec();
  double currentMass=theData->DoubleMassId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::M_PIPISCAT_NAME));
  _tMatr->evalMatrix(currentMass);

  if(_dataTypeID==1) evalElasticity(theData, currentMass);   
  else if(_dataTypeID==2){
    if(!_prodIsNotDecChannel) evalPhase(theData, currentMass);
    else evalRelativePhase(theData, currentMass);
  }
  else if(_dataTypeID==3){
    evalArgandUnits(theData, currentMass);
  }
  else{
    Alert << "_dataTypeID = " <<_dataTypeID << " is not supported!!!" << endmsg;
    exit(1); 
  }

  // else{
  //   complex<double> currentTijRel=(*_tMatr)(_prodProjectionIndex,_decProjectionIndex);
  //   complex<double> SijRel=complex<double>(1.,0.)+2.*PawianConstants::i*sqrt(thePhpVecs[_prodProjectionIndex]->factor(currentMass).real()*thePhpVecs[_decProjectionIndex]->factor(currentMass).real())*currentTijRel;
    
  //   //phase
    
  //   //note: this is a workaround
  //   if(_prodProjectionIndex!=_decProjectionIndex){
  //     SijRel=2.*PawianConstants::i*sqrt(thePhpVecs[_prodProjectionIndex]->factor(currentMass).real()*thePhpVecs[_decProjectionIndex]->factor(currentMass).real())*currentTijRel;
  //     theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::PHIFIT_PIPISCAT_NAME))=theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::PHI_PIPISCAT_NAME));
  //   } 
  //   else{
    //   complex<double> currentTijRel_rho=currentTijRel*thePhpVecs[_prodProjectionIndex]->factor(currentMass).real();
    //   double currentReERel = currentTijRel_rho.real();
    //   double currentImERel = currentTijRel_rho.imag() - 0.5;
      
    //   double phiData=theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::PHI_PIPISCAT_NAME));
    //   double deltaRel = 0.5*atan2(currentImERel, fabs(currentReERel))*PawianConstants::radToDeg + 45.0;
    //   if (currentTijRel.real()  < 0.0) {deltaRel = 180.0 - deltaRel;}
      
    //   while( (phiData-deltaRel) > 90.) deltaRel+=180.;
    //   while( (deltaRel-phiData) > 90.) deltaRel-=180.;
    //   theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::PHIFIT_PIPISCAT_NAME))=deltaRel;
    // }
    // if(thePhpVecs[_decProjectionIndex]->factor(currentMass).real()<1.e-12){//protection against calculations below threshold
    //   theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::ETAFIT_PIPISCAT_NAME))=0.;
    //   theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::PHIFIT_PIPISCAT_NAME))=0.;
    // }
    // else{
  //     double normSijRel=sqrt(norm(SijRel));
  //     theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::ETAFIT_PIPISCAT_NAME))=sqrt(norm(SijRel));
  //   }
  // }
  return (*_tMatr)(_prodProjectionIndex,_decProjectionIndex);
}

void TMatrixDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
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

  if(_dataTypeID==2 && _prodIsNotDecChannel){
    _paramNameRelPhase="deltaRel"+_kMatName+"to"+_gFactorNames.at(_decProjectionIndex);
    fitPar->Add(_paramNameRelPhase, 90., 10.);
    fitPar->SetLimits(_paramNameRelPhase, 0., 360.);
    _currentRelPhase = fitPar->Value(_paramNameRelPhase);
  }
}

void TMatrixDynamics::fillParamNameList(){
  _paramNameList.clear();

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
}

bool TMatrixDynamics::checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld){
  return true;
}

void TMatrixDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){

  //pole positions
  for(unsigned int i=0; i<_currentPoleMasses.size(); ++i){
    std::string currentPoleName=_poleNames.at(i)+"Mass";
    double currentPoleMass=fitPar->Value(currentPoleName);
    _currentPoleMasses.at(i)=currentPoleMass;
    _kPoles.at(i)->updatePoleMass(currentPoleMass);
  }
 
 //g-factors
  for(unsigned int i=0; i<_poleNames.size(); ++i){
    std::vector<double>& currentgFactorVec=_currentgFactorMap.at(i);
    for(unsigned int j=0; j<currentgFactorVec.size(); ++j){
      std::string currentName=_poleNames.at(i)+_gFactorNames.at(j)+"gFactor";
      currentgFactorVec.at(j)=fitPar->Value(currentName);
    }
    _kPoles.at(i)->updategFactors(currentgFactorVec);
  }
    
  //k-matrix bg-terms
  if(_orderKMatBg>=0){
    for(unsigned int i=0; i<=fabs(_orderKMatBg); ++i){
      for(unsigned int j=0; j<_phpVecs.size(); ++j){
	for(unsigned int k=j; k<_phpVecs.size(); ++k){
	  double newVal=fitPar->Value(_bgTermNames.at(i).at(j).at(k));
	  _currentBgTerms.at(i).at(j).at(k)=newVal;
	  _kMatr->updateBgTerms(i,j,k,newVal);
	  //	  InfoMsg << "setting current bg term(" << i << "," << j << ","<< k << ") = " << newVal <<"\tname" << _bgTermNames.at(i).at(j).at(k) << endmsg;
	}
      }
    }

    //Adler-term
    if(_withKMatAdler){
      _currentAdler0=fitPar->Value("s0"+_kMatName);
      _kMatr->updates0Adler(_currentAdler0);
    }
  }

  if(_dataTypeID==2 && _prodIsNotDecChannel){
    _currentRelPhase=fitPar->Value(_paramNameRelPhase);
  }
}



void TMatrixDynamics::init(){
  _kMatName=_kMatrixParser->keyName();
  _orderKMatBg=_kMatrixParser->orderBg();
  _withKMatAdler=_kMatrixParser->useAdler();  

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
  std::map<std::pair<std::string, std::string>, std::string> phpDescriptionVec=_kMatrixParser->phpDescriptionMap();

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

  std::pair<std::string, std::string> currentParticlePair=make_pair(firstParticleName, secondParticleName);
    std::string currentPhpDescription = phpDescriptionVec.at(currentParticlePair);

    std::vector<double> fsParticleMasses;
    fsParticleMasses.push_back(firstParticle->mass());
    fsParticleMasses.push_back(secondParticle->mass()); 

    std::shared_ptr<AbsPhaseSpace> currentPhp=PhaseSpaceFactory::instance()->getPhpPointer(currentPhpDescription, fsParticleMasses);
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
      _currentgFactorMap[i].push_back(currentGValue);
    }
  }


  std::map<int, std::vector<double> >::iterator itgFac;
  for (itgFac=_currentgFactorMap.begin(); itgFac!=_currentgFactorMap.end(); ++itgFac){
    std::vector<double> currentgVector=itgFac->second;
    std::shared_ptr<KPole> currentPole;
    if (_kMatrixParser->useBarrierFactors()) currentPole=std::shared_ptr<KPole>(new KPoleBarrier(currentgVector, _currentPoleMasses.at(itgFac->first), _phpVecs, _kMatrixParser->orbitalMom(), _kMatrixParser->useTruncatedBarrierFactors()));
    else currentPole=std::shared_ptr<KPole>(new KPole(currentgVector, _currentPoleMasses.at(itgFac->first)));
    _kPoles.push_back(currentPole);
  }

  if(_orderKMatBg<0) _kMatr=std::shared_ptr<KMatrixRel>(new KMatrixRel(_kPoles,_phpVecs ));
  else {
    _currentBgTerms.resize(_orderKMatBg+1);
    _bgTermNames.resize(_orderKMatBg+1);
    for(unsigned int i=0; i<= fabs(_orderKMatBg); ++i){
      _currentBgTerms.at(i).resize(_phpVecs.size());
      _bgTermNames.at(i).resize(_phpVecs.size());
      for(unsigned int j=0; j<_phpVecs.size(); ++j){
	_currentBgTerms.at(i).at(j).resize(_phpVecs.size());
	_bgTermNames.at(i).at(j).resize(_phpVecs.size());
	for(unsigned int k=j; k<_phpVecs.size(); ++k){
	  _currentBgTerms.at(i).at(j).at(k)=0.;

	  std::stringstream keyOrderStrStr;
	  keyOrderStrStr << i << j << k;
	  std::string keyOrder=keyOrderStrStr.str();
	  std::string currentName="bg"+keyOrder+_kMatName;
	  _bgTermNames.at(i).at(j).at(k)=currentName;
	}
      }
    }

    _kMatr=std::shared_ptr<KMatrixRel>(new KMatrixRelBg(_kPoles,_phpVecs, _orderKMatBg, _withKMatAdler ));
    if(_withKMatAdler){
      _currentAdler0=_kMatrixParser->s0Adler();
      _kMatr->updates0Adler(_currentAdler0);
      _kMatr->updatesnormAdler(_kMatrixParser->snormAdler());
    }
  }
 
  const std::string porjectionParticleNames=_kMatrixParser->projection();
  std::istringstream projParticles(porjectionParticleNames);
  std::string firstProjParticleName;
  std::string secondProjParticleName;
  projParticles >> firstProjParticleName >> secondProjParticleName;
  std::string projKey=firstProjParticleName+secondProjParticleName;    

  bool found=false;
  for(unsigned int i=0; i<_gFactorNames.size();++i){
    if(projKey==_gFactorNames[i]){
      _decProjectionIndex=i;
      found=true;
    }
  }
  if (!found){
    Alert << "decay projection index for key " << projKey << " not found" << endmsg;
    exit(0);
  }

  _tMatr=std::shared_ptr<TMatrixRel>(new TMatrixRel(_kMatr));

}

void TMatrixDynamics::evalElasticity(EvtData* theData, double currentMass){
  vector<std::shared_ptr<AbsPhaseSpace> > thePhpVecs=_tMatr->kMatrix()->phaseSpaceVec();
  complex<double> currentTijRel=(*_tMatr)(_prodProjectionIndex,_decProjectionIndex);
  complex<double> SijRel;

  //note: this is a workaround
  if(_prodProjectionIndex == _decProjectionIndex){
    SijRel=complex<double>(1.,0.)+2.*PawianConstants::i*sqrt(thePhpVecs[_prodProjectionIndex]->factor(currentMass).real()*thePhpVecs[_decProjectionIndex]->factor(currentMass).real())*currentTijRel;
  }
  else{
    SijRel=2.*PawianConstants::i*sqrt(thePhpVecs[_prodProjectionIndex]->factor(currentMass).real()*thePhpVecs[_decProjectionIndex]->factor(currentMass).real())*currentTijRel;
    }
  theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::FIT_PIPISCAT_NAME))=sqrt(norm(SijRel));
}


void TMatrixDynamics::evalPhase(EvtData* theData, double currentMass){
  vector<std::shared_ptr<AbsPhaseSpace> > thePhpVecs=_tMatr->kMatrix()->phaseSpaceVec();
  complex<double> currentTijRel=(*_tMatr)(_prodProjectionIndex,_decProjectionIndex);

  complex<double> currentTijRel_rho=currentTijRel*thePhpVecs[_prodProjectionIndex]->factor(currentMass).real();
  double currentReERel = currentTijRel_rho.real();
  double currentImERel = currentTijRel_rho.imag() - 0.5;

  double deltaRel = 0.5*atan2(currentImERel, fabs(currentReERel))*PawianConstants::radToDeg + 45.0;
  if (currentTijRel.real()  < 0.0) {deltaRel = 180.0 - deltaRel;}

  double phiData=theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::DATA_PIPISCAT_NAME));

  while( (phiData-deltaRel) > 90.) deltaRel+=180.;
  while( (deltaRel-phiData) > 90.) deltaRel-=180.;
  theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::FIT_PIPISCAT_NAME))=deltaRel;

}

void TMatrixDynamics::evalRelativePhase(EvtData* theData, double currentMass){
  vector<std::shared_ptr<AbsPhaseSpace> > thePhpVecs=_tMatr->kMatrix()->phaseSpaceVec();
  complex<double> currentTiiRel=(*_tMatr)(_prodProjectionIndex, _prodProjectionIndex);
  complex<double> currentTjjRel=(*_tMatr)(_decProjectionIndex,_decProjectionIndex);

  complex<double> currentTiiRel_rho= currentTiiRel*thePhpVecs[_prodProjectionIndex]->factor(currentMass).real();
  complex<double> currentTjjRel_rho= currentTjjRel*thePhpVecs[_decProjectionIndex]->factor(currentMass).real();

  double currentReEiiRel = currentTiiRel_rho.real();
  double currentImEiiRel = currentTiiRel_rho.imag() - 0.5;
  double deltaiiRel = 0.5*atan2(currentImEiiRel, fabs(currentReEiiRel))*PawianConstants::radToDeg + 45.0;
  if (currentReEiiRel  < 0.0) {deltaiiRel = 180.0 - deltaiiRel;}

  double currentReEjjRel = currentTjjRel_rho.real();
  double currentImEjjRel = currentTjjRel_rho.imag() - 0.5;
  double deltajjRel = 0.5*atan2(currentImEjjRel, fabs(currentReEjjRel))*PawianConstants::radToDeg + 45.0;
  if (currentReEjjRel  < 0.0) {deltajjRel = 180.0 - deltajjRel;}

  theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::FIT_PIPISCAT_NAME))=deltaiiRel+deltajjRel+_currentRelPhase;
}

void TMatrixDynamics::evalArgandUnits(EvtData* theData, double currentMass){
  vector<std::shared_ptr<AbsPhaseSpace> > thePhpVecs=_tMatr->kMatrix()->phaseSpaceVec();
  complex<double> currentTijRel=(*_tMatr)(_prodProjectionIndex,_decProjectionIndex);

  double sqrTij=0.;
  if( thePhpVecs[_prodProjectionIndex]->factor(currentMass).real() > 1.e-10 && thePhpVecs[_decProjectionIndex]->factor(currentMass).real() > 1.e-10){
    sqrTij=norm(sqrt(thePhpVecs[_prodProjectionIndex]->factor(currentMass).real())*currentTijRel*sqrt(thePhpVecs[_decProjectionIndex]->factor(currentMass).real()));
  }

  theData->DoubleId.at(IdStringMapRegistry::instance()->stringId(EvtDataScatteringList::FIT_PIPISCAT_NAME))=sqrTij;
}


