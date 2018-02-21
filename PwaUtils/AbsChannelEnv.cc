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

// AbsChannelEnv class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy


#include "AbsChannelEnv.hh"
#include "ConfigParser/ParserBase.hh"
#include "GlobalEnv.hh"
#include "AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/IsobarTensorDecay.hh"
#include "PwaUtils/IsobarTensorPsiToGamXDecay.hh"
#include "PwaUtils/OmegaTo3PiLSDecay.hh"
#include "PwaUtils/OmegaTo3PiTensorDecay.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/ProdChannelInfo.hh"
#include "PwaUtils/FsParticleProjections.hh"
#include "PwaUtils/AbsHist.hh"
#include "PwaUtils/RootHist.hh"
#include "PwaUtils/PhpGenDynamics.hh"
#include "Particle/Particle.hh"
#include "Particle/PdtParser.hh"
#include "FitParams/ParamFactory.hh" 
#include "Event/MassRangeCut.hh"
#include "ErrLogger/ErrLogger.hh"



short AbsChannelEnv::CHANNEL_PBARP = 1;
short AbsChannelEnv::CHANNEL_EPEM = 2;
short AbsChannelEnv::CHANNEL_RES = 3;
short AbsChannelEnv::CHANNEL_GAMMAP = 4;
short AbsChannelEnv::CHANNEL_PIPISCATTERING = 5;
short AbsChannelEnv::CHANNEL_GG = 6;

AbsChannelEnv::AbsChannelEnv(ParserBase* theParser, short channelType) :
  _channelType(channelType)
  ,_alreadySetUp(false)
  ,_noFinalStateParticles(0)
  ,_absDecList(new AbsDecayList())
  ,_prodDecList(new AbsDecayList())
  ,_useMassRange(false)
  ,_usePhpGenDynamics(false)
  ,_theParser(theParser)
  ,_cmEnergy(0.)
  ,_channelSubTypeName("")
  ,_initial4Vec(Vector4<double>(0.,0.,0.,0.))
  ,_genRange(std::pair<float,float>(0.,0.))
  ,_projectile4Vec(Vector4<double>(0.,0.,0.,0.))
{
  std::string init4VecStr=theParser->intitial4Vec();
  _initial4Vec=extract4Vec(init4VecStr);
  InfoMsg << "initial4Vec: " << _initial4Vec << endmsg;

  std::string projectile4VecStr=theParser->projectile4Vec();
  _projectile4Vec=extract4Vec(projectile4VecStr);
  InfoMsg << "projectile4Vec: " << _projectile4Vec << endmsg;
  _theLh.reset();
}



std::shared_ptr<AbsLh> AbsChannelEnv::Lh(){
   if(_theLh == nullptr){
      Alert << "Accessing Lh NULL pointer." << endmsg;
   }
   return _theLh;
}



void AbsChannelEnv::setupGlobal(ChannelID id){
   if(_alreadySetUp){
      Alert << "ChannelEnv already set up!" << endmsg;
      exit(1);
   }
   _alreadySetUp = true;

   // Set channel id
   _channelID = id;

   // Event weights for data?
   _useDataEvtWeight = _theParser->useDataEvtWeight();
   _useMCEvtWeight = _theParser->useMCEvtWeight();

   // cloned particles
   const std::vector<std::string> cloneParticle=_theParser->cloneParticle();
   std::vector<std::string>::const_iterator itcP;

   for ( itcP = cloneParticle.begin(); itcP != cloneParticle.end(); ++itcP){
      std::istringstream particles(*itcP);
      std::string particleOld;
      std::string particleNew;
      particles >> particleOld >> particleNew;

      GlobalEnv::instance()->particleTable()->clone(particleNew, particleOld);
   }


   //final state particles
   const std::vector<std::string> finalStateParticleStr=_theParser->finalStateParticles();

   std::vector<std::string>::const_iterator itStr;
   for ( itStr = finalStateParticleStr.begin(); itStr != finalStateParticleStr.end(); ++itStr){
      Particle* currentParticle = GlobalEnv::instance()->particleTable()->particle(*itStr);
      _finalStateParticles.push_back(currentParticle);
   }

   _noFinalStateParticles= (int) _finalStateParticles.size();


   _fsParticleProjections = std::shared_ptr<FsParticleProjections>(new FsParticleProjections(_finalStateParticles)); 
  //decays

  std::vector<std::string> decaySystem= _theParser->decaySystem();
  for ( itStr = decaySystem.begin(); itStr != decaySystem.end(); ++itStr){

    Particle* motherParticle =0;
    std::vector<Particle*> daughterParticles;

    std::stringstream stringStr;
    stringStr << (*itStr);

    std::string tmpName;

    bool useIsospin=true;
    bool isDecParticle=false;
    bool firstArgument=true;
    bool secondArgument=false;
    std::string usedSystem("");

    while(stringStr >> tmpName){
      if(firstArgument){
	usedSystem=tmpName;
	firstArgument=false;
	secondArgument=true;
	continue;
      }
      if(secondArgument){
	if(tmpName=="noIso"){
	 useIsospin=false;
         continue;
        }
	secondArgument=false;
      }
      if(tmpName=="To") {
        isDecParticle=true;
        continue;
      }
      if (GlobalEnv::instance()->particleTable()->particle(tmpName) == NULL){
	Alert << "Particle " << tmpName << " does not exist in ParticleTable. Please clone an existing particle or add it to the table.";
	exit(1);
      }
      if(isDecParticle){
	daughterParticles.push_back(GlobalEnv::instance()->particleTable()->particle(tmpName));
      }
      else{
	motherParticle = GlobalEnv::instance()->particleTable()->particle(tmpName);
      }
    }
    std::shared_ptr<AbsDecay> tmpDec;
    if(daughterParticles.size()==2){
      if (usedSystem=="Heli") tmpDec= std::shared_ptr<AbsDecay>(new IsobarHeliDecay(motherParticle, daughterParticles[0], daughterParticles[1], _channelID));
      else if (usedSystem=="Cano")  tmpDec= std::shared_ptr<AbsDecay>(new IsobarLSDecay(motherParticle, daughterParticles[0], daughterParticles[1], _channelID));
      else if (usedSystem=="Tensor"){
	if(motherParticle->twoJ()==2 && motherParticle->theParity()==-1 && motherParticle->theCParity()==-1){
	  if (daughterParticles[0]->name()=="photon") tmpDec= std::shared_ptr<AbsDecay>(new IsobarTensorPsiToGamXDecay(motherParticle, daughterParticles[0], daughterParticles[1], _channelID));
	  else if (daughterParticles[1]->name()=="photon") tmpDec= std::shared_ptr<AbsDecay>(new IsobarTensorPsiToGamXDecay(motherParticle, daughterParticles[1], daughterParticles[0], _channelID));
	  else tmpDec= std::shared_ptr<AbsDecay>(new IsobarTensorDecay(motherParticle, daughterParticles[0], daughterParticles[1], _channelID));
	}
	else tmpDec= std::shared_ptr<AbsDecay>(new IsobarTensorDecay(motherParticle, daughterParticles[0], daughterParticles[1], _channelID));
      }
      else {
	Alert << "used decay system\t" << usedSystem << "\tnot supported!!!\n" << endmsg;
	exit(1);
      }

    }

    else if(daughterParticles.size()==3){
      if (usedSystem=="Cano") tmpDec= std::shared_ptr<AbsDecay>(new OmegaTo3PiLSDecay(motherParticle, daughterParticles[0], daughterParticles[1], daughterParticles[2], _channelID));
      else if (usedSystem=="Tensor") tmpDec= std::shared_ptr<AbsDecay>(new OmegaTo3PiTensorDecay(motherParticle, daughterParticles[0], daughterParticles[1], daughterParticles[2], _channelID));
      else {
	Alert << "used decay system\t" << usedSystem << "\tnot supported!!!\n" << endmsg;
	exit(1);
      }
    }

    else {
      Alert << "Decay\t" << (*itStr) << "\tnot supported!!!" ;
    }

    if(!useIsospin) tmpDec->disableIsospin();
    tmpDec->extractStates();
    _absDecList->addDecay(tmpDec);
    InfoMsg << "added decay " << tmpDec->name() << " to decay list" << endmsg;
  }


  //produce particle pairs
  std::vector<std::string> productionSystem = _theParser->productionSystem();

  for ( itStr = productionSystem.begin(); itStr != productionSystem.end(); ++itStr){
    std::string currentLine= (*itStr);
    std::shared_ptr<ProdChannelInfo> currentChannelInfo(new ProdChannelInfo(currentLine) );
    currentChannelInfo->print(std::cout);
    _prodChannelInfoList.push_back(currentChannelInfo);
  }

  //set prefactor
  std::vector<std::string> preFactorVec = _theParser->preFactor();
  for ( itStr = preFactorVec.begin(); itStr != preFactorVec.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);

    std::string ampName;
    stringStr >> ampName;

    std::string valueStr;
    stringStr >> valueStr;

    double currentValue;
    currentValue=atof(valueStr.c_str());

    _preFactorMap[ampName]=currentValue;
  }


  //fill vector histMassSystems
  std::vector<std::string> theHistMassNames=_theParser->histMassNames();
  for ( itStr = theHistMassNames.begin(); itStr != theHistMassNames.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);

    std::string tmpName;
    std::vector<std::string> currentStringVec;
    while(stringStr >> tmpName){
      currentStringVec.push_back(tmpName);
    }
    _histMassSystems.push_back(currentStringVec);
  }

  //mass range cuts
  std::vector<std::string> massRangeCuts=_theParser->massRangeCuts();
  if(massRangeCuts.size()>0) _useMassRange=true;
  for ( itStr = massRangeCuts.begin(); itStr != massRangeCuts.end(); ++itStr){
    std::string currentString=*itStr;
    std::shared_ptr<MassRangeCut> currentMassRangeCut(new MassRangeCut( currentString, _finalStateParticles));
    _massRangeCuts.push_back(currentMassRangeCut);
  }

  //generator range
  std::string genRange=_theParser->genRange();
  std::string _minGenMass, _maxGenMass;
  std::istringstream iss(genRange);
  iss >> _minGenMass;
  iss >> _maxGenMass;
  _genRange = std::make_pair(atof(_minGenMass.c_str()), atof(_maxGenMass.c_str()));

  //php gen dynamics
  std::vector<std::string> phpGenDynamics=_theParser->phpGenDynamics();
  if(phpGenDynamics.size()>0) _usePhpGenDynamics=true;
  for ( itStr = phpGenDynamics.begin(); itStr != phpGenDynamics.end(); ++itStr){
    std::string currentString=*itStr;
    std::shared_ptr<PhpGenDynamics> currentPhpGenDynamic(new PhpGenDynamics( currentString, _finalStateParticles));
    _phpGenDynamics.push_back(currentPhpGenDynamic);
  }  

 // hist angles
  std::vector<std::string> theHistAngleNames=_theParser->histAngleNames();
//fill vector histMassSystems
  for ( itStr = theHistAngleNames.begin(); itStr != theHistAngleNames.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);

    std::string tmpName;
    std::vector<std::string> currentStringDecVec;
    std::vector<std::string> currentStringDecVec2;
    std::vector<std::string> currentStringMotherVec;
    bool isDecParticle=true;
    bool isFirstDecParticle=true;
    short nBodyDecay=2;
    while(stringStr >> tmpName){
      if(tmpName=="from") {
	isDecParticle=false;
	continue;
      }
      else if(tmpName=="and") {
	isFirstDecParticle=false;
	nBodyDecay=3;
	continue;
      }
      if(isDecParticle && isFirstDecParticle) currentStringDecVec.push_back(tmpName);
      else if(isDecParticle && !isFirstDecParticle) currentStringDecVec2.push_back(tmpName);
      else currentStringMotherVec.push_back(tmpName);
    }
    std::shared_ptr<angleHistData> currentAngleHistData(new angleHistData(currentStringMotherVec, currentStringDecVec, currentStringDecVec2, nBodyDecay));
    _angleHistDataVec.push_back(currentAngleHistData);
  }

  // 2D hists
  std::vector<std::string> theHistAngleNames2D=_theParser->histAngleNames2D();
  for ( itStr = theHistAngleNames2D.begin(); itStr != theHistAngleNames2D.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);

    std::string tmpName;
    std::vector<std::string> currentStringDecVec;
    std::vector<std::string> currentStringMotherVec;
    std::vector<std::string> currentStringDecVec_2;
    std::vector<std::string> currentStringMotherVec_2;
    bool isDecParticle=true;
    bool isFirstAxes=true;
    while(stringStr >> tmpName){
      if(tmpName=="from") {
        isDecParticle=false;
	std::cout << tmpName << std::endl;
        continue;
      }
      else if(tmpName=="versus") {
        isFirstAxes=false;
	isDecParticle=true;
        continue;
      }
      if(isDecParticle && isFirstAxes) currentStringDecVec.push_back(tmpName);
      else if (!isDecParticle && isFirstAxes) currentStringMotherVec.push_back(tmpName);
      else if(isDecParticle && !isFirstAxes) currentStringDecVec_2.push_back(tmpName);
      else currentStringMotherVec_2.push_back(tmpName);
    }
    std::shared_ptr<angleHistData2D> currentAngleHistData2D(new angleHistData2D(currentStringMotherVec, currentStringDecVec, currentStringMotherVec_2, currentStringDecVec_2));
    _angleHistDataVec2D.push_back(currentAngleHistData2D);
  }

 // calculate contributions
  std::vector<std::string> theCalcContribution=_theParser->calcContribution();
  for ( itStr = theCalcContribution.begin(); itStr != theCalcContribution.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);

    std::string tmpName;
    std::string contribName;
    std::vector<std::string> currentStringZeroAmp;
    bool isContribName = true;
    while(stringStr >> tmpName){
      if(tmpName=="withZeroAmp") {
      	isContribName = false;
      	continue;
      }
      if(isContribName) contribName = tmpName;
      else currentStringZeroAmp.push_back(tmpName);
    }
    std::shared_ptr<calcContributionData> currentCalcContributionData(new calcContributionData(contribName, currentStringZeroAmp));
    _calcContributionDataVec.push_back(currentCalcContributionData);
  }

  GlobalEnv::instance()->fillReplacedSuffixMap(_theParser);
}

std::shared_ptr<AbsPawianParameters> AbsChannelEnv::defaultPawianParams(){
  std::shared_ptr<AbsPawianParameters> result=ParamFactory::instance()->getParametersPointer("Pawian");
  Lh()->fillDefaultParams(result);
  return result;
}

void AbsChannelEnv::setWignerDRefs(){
  std::vector<std::shared_ptr<AbsDecay> > theProdVec=_prodDecList->getList();
  std::vector<std::shared_ptr<AbsDecay> >::iterator it;
  std::string prodKey="prod";
  for (it=theProdVec.begin(); it!=theProdVec.end(); ++it){
    (*it)->setWigDRefKey(prodKey);
  }
}

bool AbsChannelEnv::checkReactionChain(){
  InfoMsg << "check complete reaction chains" << endmsg;

  // first check: unique names for final state particles
  for(unsigned int i=0; i<_finalStateParticles.size()-1; ++i){
    Particle* currentParticlei=_finalStateParticles.at(i);
    for(unsigned int j=i+1; j<_finalStateParticles.size(); ++j){
      Particle* currentParticlej=_finalStateParticles.at(j);
      if(currentParticlei->name() == currentParticlej->name()){
	Alert << "Unique names are requested for particles in the final state particle list!!!!"
	      << "\n there are at least two particles with the name: " << currentParticlei->name() 
	      << endmsg;
	exit(1);
      }
    }
  }

  // second check: each decay amplitude must contain at least one mother amplitude and at least on production amplitude
  std::vector<std::shared_ptr<AbsDecay> >::iterator it;
  std::vector<std::shared_ptr<AbsDecay> > theDecVec=_absDecList->getList();
  for (it=theDecVec.begin(); it!=theDecVec.end(); ++it){
    std::vector< std::shared_ptr<AbsDecay> > prodAmpRefList=(*it)->prodAmpRefList();
    if( prodAmpRefList.size() < 1 ){
      Alert << "decay amplitude: " << (*it)->name() << " does not contain at least one production amplitude!!!!" << endmsg;
      exit(1);
    }
    std::vector< std::shared_ptr<AbsDecay> > motherAmpRefList=(*it)->motherAmpRefList();
    if( motherAmpRefList.size() < 1 ){
      Alert << "decay amplitude: " << (*it)->name() << " does not contain at least one mother amplitude!!!!" << endmsg;
      exit(1);
    }
  }

  std::vector<Particle*>::iterator itParticle;
  for (itParticle=_finalStateParticles.begin(); itParticle!=_finalStateParticles.end(); ++itParticle){
    InfoMsg << "paricle name: " << (*itParticle)->name() << endmsg;
  }
 
  // third check: each production amplitude must finally decay in all final state particles
  std::vector<std::shared_ptr<AbsDecay> > theProdVec=_prodDecList->getList();
  for (it=theProdVec.begin(); it!=theProdVec.end(); ++it){
    std::vector<std::string> fsPartListNames = (*it)->fsParticleNameList();
    for (itParticle=_finalStateParticles.begin(); itParticle!=_finalStateParticles.end(); ++itParticle){
      std::vector<std::string>::const_iterator itToFind = std::find(fsPartListNames.begin(), fsPartListNames.end(), (*itParticle)->name());
      if (itToFind != fsPartListNames.end()){
	//found
       fsPartListNames.erase(itToFind);
      }
      else{
        Alert << "production amplitude " << (*it)->name() << " does finally not decay in exactly all final state particles!!!"
	      << "\nparticle with the name " << (*itParticle)->name() << " not found!!!" << endmsg;
      	exit(1);
      }
    }
  }

  InfoMsg << "reaction chains have been successfully checked!!!" << endmsg;
  return true;
}


std::shared_ptr<AbsHist> AbsChannelEnv::CreateHistInstance(std::string additionalSuffix, bool withTruth){

  return std::shared_ptr<AbsHist>(new RootHist(additionalSuffix, withTruth));
}

void AbsChannelEnv::replaceParameterSuffixes(){
  std::map<std::string, std::string> decSuffixNames=GlobalEnv::instance()->toBeReplacedSuffixMap();
  std::map<std::string, std::string>::iterator itMapStrStr;
  for (itMapStrStr=decSuffixNames.begin(); itMapStrStr!=decSuffixNames.end(); ++itMapStrStr){
    _absDecList->replaceSuffix(itMapStrStr->first, itMapStrStr->second);
    _prodDecList->replaceSuffix(itMapStrStr->first, itMapStrStr->second);
  }
}

void AbsChannelEnv::replaceMassKeys(){
  std::vector<std::string> replMassKeyVec = _theParser->replaceMassKey();
  std::map<std::string, std::string> decRepMassKeyNames;

  std::vector<std::string>::const_iterator itStr;
  for ( itStr = replMassKeyVec.begin(); itStr != replMassKeyVec.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);
    std::string oldStr;
    stringStr >> oldStr;

    std::string newStr;
    stringStr >> newStr;
    decRepMassKeyNames[oldStr]=newStr;
  }

  std::map<std::string, std::string>::iterator itMapStrStr;
  for (itMapStrStr=decRepMassKeyNames.begin(); itMapStrStr!=decRepMassKeyNames.end(); ++itMapStrStr){
    _absDecList->replaceMassKey(itMapStrStr->first, itMapStrStr->second);
  }
}

void AbsChannelEnv::addDynamics(){
  std::vector<std::shared_ptr<AbsDecay> > absDecList= _absDecList->getList();
  std::vector<std::string> decDynVec = _theParser->decayDynamics();

  std::vector<std::string>::const_iterator itStr;
  for ( itStr = decDynVec.begin(); itStr != decDynVec.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);

    std::string particleStr;
    stringStr >> particleStr;

    std::string dynStr;
    stringStr >> dynStr;

    std::string tmpName;
    std::vector<std::string> additionalStringVec;
    while(stringStr >> tmpName){
      additionalStringVec.push_back(tmpName);
    }

    std::vector<std::shared_ptr<AbsDecay> >::iterator itDec;
    for (itDec=absDecList.begin(); itDec!=absDecList.end(); ++itDec){
      std::string theDecName=(*itDec)->name();
      std::string toFind=particleStr+"To";
      size_t found;
      found = theDecName.find(toFind);

      if (found!=string::npos && found==0){
	(*itDec)->enableDynamics(dynStr, additionalStringVec);
      }
    }
  }
}

void AbsChannelEnv::setDecayLevels(){
   std::vector<std::shared_ptr<AbsDecay> > prodDecList= _prodDecList->getList();
   std::vector<std::shared_ptr<AbsDecay> >::iterator itProdDecList;
   for (itProdDecList=prodDecList.begin(); itProdDecList!=prodDecList.end(); ++itProdDecList){
     (*itProdDecList)->setDecayLevelTree(AbsDecay::decayLevel::isProdAmp, *itProdDecList, *itProdDecList);    
   }
}

void AbsChannelEnv::setPrefactors(){
 //set prefactor for production and decay amplitudes
  std::map<std::string, double>::iterator strDoubleIt;
  for(strDoubleIt=_preFactorMap.begin(); strDoubleIt!=_preFactorMap.end(); ++strDoubleIt){
    std::string currentAmplitudeName=strDoubleIt->first;
    double currentPrefactor=strDoubleIt->second;
    
    std::shared_ptr<AbsDecay> currentDec=_prodDecList->decay(currentAmplitudeName);
    if(0!=currentDec){
      currentDec->setPreFactor(currentPrefactor);
      // InfoMsg << "Set prefactor " << currentPrefactor << " for amplitude " << currentAmplitudeName << endmsg;
      currentDec->disableIsospin();
      InfoMsg << "Disable isospin coupling and set prefactor " << currentPrefactor << " for amplitude " << currentAmplitudeName << endmsg;
    }
    else{
      // look in decay amplitudes
      currentDec=_absDecList->decay(currentAmplitudeName);
      if(0!=currentDec){
	currentDec->setPreFactor(currentPrefactor);
	InfoMsg << "Set prefactor " << currentPrefactor << " for amplitude " << currentAmplitudeName << endmsg;
      }
      else{
	Alert << "Amplitude with name\t" << currentAmplitudeName << "\tnot found!!!" << endmsg;
	exit(0);
      }
    }
  }
}

Vector4<double> AbsChannelEnv::extract4Vec(std::string theString){
   std::stringstream stringStr;
  stringStr << theString;
 
  std::string eValStr;
  stringStr >> eValStr;
  double eVal;
  eVal=atof(eValStr.c_str());
 
  std::string pxValStr;
  stringStr >> pxValStr;
  double pxVal;
  pxVal=atof(pxValStr.c_str()); 

  std::string pyValStr;
  stringStr >> pyValStr;
  double pyVal;
  pyVal=atof(pyValStr.c_str());

  std::string pzValStr;
  stringStr >> pzValStr;
  double pzVal;
  pzVal=atof(pzValStr.c_str());

  Vector4<double> result(eVal, pxVal, pyVal, pzVal);
  return result;   
}

