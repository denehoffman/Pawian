//************************************************************************//
//									  //
//  Copyright 2016 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// RootHist class definition file. -*- C++ -*-
// Copyright 2016 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>

#include "PwaUtils/RootHist.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"

#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/IdStringMapRegistry.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/GlobalEnv.hh"
//#include "PwaUtils/EvtDataBaseList.hh"
#include "FitParams/AbsPawianParameters.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"

#include "TLorentzVector.h"
#include "ErrLogger/ErrLogger.hh"
#include "TTree.h"

RootHist::RootHist(std::string additionalSuffix, bool withTruth) :
  AbsHist(additionalSuffix, withTruth)
 {
  std::ostringstream rootFileName;
  rootFileName << "./pawianHists" << GlobalEnv::instance()->outputFileNameSuffix() << _additionalSuffix.c_str() <<  ".root";
  _theTFile=new TFile(rootFileName.str().c_str(),"recreate");
  
  _dataFourvecs = new TTree("_dataFourvecs", "_dataFourvecs");
  _fittedFourvecs = new TTree("_fittedFourvecs", "_fittedFourvecs");
  if (_fillTruths) _truthFourvecs = new TTree("_truthFourvecs", "_truthFourvecs");
 
  std::vector<std::shared_ptr<angleHistData> >::iterator itAngleVec;
  for (itAngleVec=_angleHistDataVec.begin(); itAngleVec!=_angleHistDataVec.end(); ++itAngleVec){

    //Heli data
    initAngleHists(_angleDataHistMap, (*itAngleVec), "Data", "Heli");

    //Gottfried Jackson data
    initAngleHists(_angleGJDataHistMap, (*itAngleVec), "Data", "GJ");

    //helicity MC
    initAngleHists(_angleMcHistMap, (*itAngleVec), "Mc", "Heli");

    //Gottfried Jackson MC
    initAngleHists(_angleGJMcHistMap, (*itAngleVec), "Mc", "GJ");

    //helicity fit
    initAngleHists(_angleFitHistMap, (*itAngleVec), "Fit", "Heli");

    //Gottfried Jackson fit
    initAngleHists(_angleGJFitHistMap, (*itAngleVec), "Fit", "GJ");

    if(_fillTruths){
    //helicity truth
    initAngleHists(_angleTruthHistMap, (*itAngleVec), "TruthWoWeight", "Heli");

    //Gottfried Jackson truth
    initAngleHists(_angleGJTruthHistMap, (*itAngleVec), "TruthWoWeight", "GJ");

    //helicity fit
    initAngleHists(_angleTruthFitHistMap, (*itAngleVec), "TruthWWeight", "Heli");

    //Gottfried Jackson fit
    initAngleHists(_angleGJTruthFitHistMap, (*itAngleVec), "TruthWWeigth", "GJ");
    }
  }

  // 2D-Histograms
  std::vector<std::shared_ptr<angleHistData2D> >::iterator itAngleVec2D;
  for (itAngleVec2D=_angleHistDataVec2D.begin(); itAngleVec2D!=_angleHistDataVec2D.end(); ++itAngleVec2D){
    std::string tmpBaseName= (*itAngleVec2D)->_name;
    boost::replace_all(tmpBaseName,"+","p");
    boost::replace_all(tmpBaseName,"-","m");
    std::string histThetaName="DataTheta"+tmpBaseName;
    std::string histPhiName="DataPhi"+tmpBaseName;
    std::string histThetaDescription = "cos#Theta(" +(*itAngleVec2D)->_name + ") (data); " + (*itAngleVec2D)->_nameXAxis + " ; " + (*itAngleVec2D)->_nameYAxis;
    std::string histPhiDescription = "#phi(" +(*itAngleVec2D)->_name + ") (data); " + (*itAngleVec2D)->_nameXAxis + " ; " + (*itAngleVec2D)->_nameYAxis;

    TH2F* currentThetaAngleDataHist=new TH2F(histThetaName.c_str(), histThetaDescription.c_str(), 100., -1., 1., 100., -1., 1.);
    TH2F* currentPhiAngleDataHist=new TH2F(histPhiName.c_str(), histPhiDescription.c_str(), 100., -3.14159, 3.14159, 100., -3.14159, 3.1415);
    currentThetaAngleDataHist->Sumw2();
    currentPhiAngleDataHist->Sumw2();
    _angleDataHistMap2D[*itAngleVec2D].push_back(currentThetaAngleDataHist);
    _angleDataHistMap2D[*itAngleVec2D].push_back(currentPhiAngleDataHist);

    histThetaName="MCTheta"+tmpBaseName;
    histPhiName="MCPhi"+tmpBaseName;
    histThetaDescription = "cos#Theta(" +(*itAngleVec2D)->_name + ") (MC); " + (*itAngleVec2D)->_nameXAxis + " ; " + (*itAngleVec2D)->_nameYAxis;
    histPhiDescription = "#phi(" +(*itAngleVec2D)->_name + ") (MC); " + (*itAngleVec2D)->_nameXAxis + " ; " + (*itAngleVec2D)->_nameYAxis;

    TH2F* currentThetaAngleMcHist=new TH2F(histThetaName.c_str(), histThetaDescription.c_str(), 100., -1., 1., 100., -1., 1.);
    TH2F* currentPhiAngleMcHist=new TH2F(histPhiName.c_str(), histPhiDescription.c_str(), 100., -3.14159, 3.14159, 100., -3.14159, 3.14159);
    currentThetaAngleMcHist->Sumw2();
    currentPhiAngleMcHist->Sumw2();
    _angleMcHistMap2D[*itAngleVec2D].push_back(currentThetaAngleMcHist);
    _angleMcHistMap2D[*itAngleVec2D].push_back(currentPhiAngleMcHist);

    histThetaName="FitTheta"+tmpBaseName;
    histPhiName="FitPhi"+tmpBaseName;
    histThetaDescription = "cos#Theta(" +(*itAngleVec2D)->_name + ") (fit); " + (*itAngleVec2D)->_nameXAxis + " ; " + (*itAngleVec2D)->_nameYAxis;
    histPhiDescription = "#phi(" +(*itAngleVec2D)->_name + ") (fit); " + (*itAngleVec2D)->_nameXAxis + " ; " + (*itAngleVec2D)->_nameYAxis;

    TH2F* currentThetaAngleFitHist=new TH2F(histThetaName.c_str(), histThetaDescription.c_str(), 100., -1., 1., 100., -1., 1.);
    TH2F* currentPhiAngleFitHist=new TH2F(histPhiName.c_str(), histPhiDescription.c_str(), 100., -3.14159, 3.14159, 100., -3.14159, 3.14159);
    currentThetaAngleFitHist->Sumw2();
    currentPhiAngleFitHist->Sumw2();
    _angleFitHistMap2D[*itAngleVec2D].push_back(currentThetaAngleFitHist);
    _angleFitHistMap2D[*itAngleVec2D].push_back(currentPhiAngleFitHist);
  }


  //invariant masses
  std::vector<std::vector<std::string> > histMassNameVec=GlobalEnv::instance()->Channel()->histMassSystems();
  std::vector<std::vector<std::string> >::iterator itVecStr;
  for(itVecStr=histMassNameVec.begin(); itVecStr!=histMassNameVec.end(); ++itVecStr){
    std::shared_ptr<massHistData> tmpMassHistData(new massHistData(*itVecStr));
    std::string tmpBaseName=tmpMassHistData->_name;
    boost::replace_all(tmpBaseName,"+","p");
    boost::replace_all(tmpBaseName,"-","m");
    std::string histName="Data"+tmpBaseName;
    std::string histDescription = "M("+tmpMassHistData->_name+") (data)";

        double massMin = 0.;
    double massMax=_cmMass;

      
    std::vector<std::string> fspNames=tmpMassHistData->_fspNames;
    std::vector<Particle*> allFsp = GlobalEnv::instance()->Channel()->finalStateParticles();
    std::vector<Particle*>::iterator itAllFsp;

        for(itAllFsp = allFsp.begin(); itAllFsp != allFsp.end(); ++itAllFsp){
       bool isObserver = true;
       std::vector<std::string>::iterator itStr2;
       for(itStr2=fspNames.begin(); itStr2!=fspNames.end(); ++itStr2){
	  if(*itStr2 == (*itAllFsp)->name())
	     isObserver = false;
       }
       if(isObserver)
	  massMax -= (*itAllFsp)->mass();
       else
	  massMin += (*itAllFsp)->mass();
    }

    massMax += (massMax - massMin) * 0.02;
    massMin -= (massMax - massMin) * 0.02;

        TH1F* currentMassDataHist=new TH1F(histName.c_str(), histDescription.c_str(), 100., massMin, massMax);
    currentMassDataHist->Sumw2();
    _massDataHistMap[tmpMassHistData]=currentMassDataHist;

    histName="MC"+tmpBaseName;
    histDescription = "M("+tmpMassHistData->_name+") (MC)";
    TH1F* currentMassMcHist=new TH1F(histName.c_str(), histDescription.c_str(), 100., massMin, massMax);
    currentMassMcHist->Sumw2();
    _massMcHistMap[tmpMassHistData]=currentMassMcHist;

    histName="Fit"+tmpBaseName;
    histDescription = "M("+tmpMassHistData->_name+") (fit)";
    TH1F* currentMassFitHist=new TH1F(histName.c_str(), histDescription.c_str(), 100., massMin, massMax);
    currentMassFitHist->Sumw2();
    _massFitHistMap[tmpMassHistData]=currentMassFitHist;

    if (_fillTruths){
      histName="TruthWoWeight"+tmpBaseName;
      histDescription = "M("+tmpMassHistData->_name+") (truthWoWeight)";
      TH1F* currentMassTruthHist=new TH1F(histName.c_str(), histDescription.c_str(), 100., massMin, massMax);
      currentMassTruthHist->Sumw2();
      _massTruthHistMap[tmpMassHistData]=currentMassTruthHist;
      
      histName="TruthWWeight"+tmpBaseName;
      histDescription = "M("+tmpMassHistData->_name+") (truthWWeight)";
      TH1F* currentMassTruthFitHist=new TH1F(histName.c_str(), histDescription.c_str(), 100., massMin, massMax);
      currentMassTruthFitHist->Sumw2();
      _massTruthFitHistMap[tmpMassHistData]=currentMassTruthFitHist;
    }
  }
  
  //tree
  for(auto fsIt = _fsParticles.begin(); fsIt != _fsParticles.end(); ++fsIt){
     std::string particleName = (*fsIt)->name();
     _fourVecMap[particleName] = std::shared_ptr<TLorentzVector>(new TLorentzVector(0,0,0,0));
     _dataFourvecs->Branch(particleName.c_str(), "TLorentzVector", _fourVecMap[particleName].get());
     _fittedFourvecs->Branch(particleName.c_str(), "TLorentzVector", _fourVecMap[particleName].get());
     if(_fillTruths) _truthFourvecs->Branch(particleName.c_str(), "TLorentzVector", _fourVecMap[particleName].get()); 
  }

  _dataFourvecs->Branch("weight", &_weightToWrite, "weight");
  _fittedFourvecs->Branch("weight", &_weightToWrite, "weight");
  if(_fillTruths) _truthFourvecs->Branch("weight", &_weightToWrite, "weight");
 }

RootHist::~RootHist(){
  _theTFile->Write();
  _theTFile->Close();
}

void RootHist::fillFromLhData(std::shared_ptr<AbsLh> theLh, std::shared_ptr<AbsPawianParameters> fitParams){

  if(0==theLh){
    Alert <<"AbsLh* is a 0 pointer !!!!" ;  // << endmsg;
    exit(1);
  }

  theLh->updateFitParams(fitParams);

  const std::vector<EvtData*> dataList=theLh->getDataVec();
  double integralDataWWeight=0.;

  std::vector<EvtData*>::const_iterator it=dataList.begin();
  while(it!=dataList.end())
    {
      double weight = (*it)->evtWeight;
      integralDataWWeight+=weight;
      fillEvt((*it), weight, "data");
      ++it;
    }

  //  const std::vector<EvtData*> mcList=theEvtList->getMcVecs();
  const std::vector<EvtData*> mcList=theLh->getMcVec();
  double integralMC=0.;
  //  double integralFitWeight=0.;

  it=mcList.begin();
  while(it!=mcList.end())
    {
      double evtWeight = (*it)->evtWeight;
      integralMC+=evtWeight;
      fillEvt((*it), evtWeight, "mc");

      double fitWeight= theLh->calcEvtIntensity( (*it), fitParams );
      //      integralFitWeight+=fitWeight;
      fillEvt((*it), evtWeight*fitWeight, "fit");

      ++it;
    }

  //  std::string outputFileName= "qaSummary" + GlobalEnv::instance()->outputFileNameSuffix() + ".dat";
  //  std::ofstream theQaStream ( outputFileName.c_str(), std::ofstream::out | std::ofstream::app);

  //  double integralDataWoWeight=(double) dataList.size();

  //  InfoMsg <<"No of data events without weight " << integralDataWoWeight ;  // << endmsg;
  //  theQaStream <<"No of data events without weight " << integralDataWoWeight << "\n";
  //
  //  InfoMsg <<"No of data events with weight " << integralDataWWeight ;  // << endmsg;
  //  theQaStream <<"No of data events with weight " << integralDataWWeight << "\n";
  //
  //  InfoMsg <<"No of MC events " << integralMC ;  // << endmsg;
  //  theQaStream <<"No of MC events " << integralMC << "\n";
  //
  double scaleFactor = integralDataWWeight/integralMC;
  //  InfoMsg <<"scaling factor  " << scaleFactor;  // << endmsg;
  //  theQaStream <<"scaling factor  " << scaleFactor << "\n";  // << endmsg;
  //
  //  InfoMsg <<"no of fitted events with scaling factor: " << integralFitWeight*scaleFactor ;  // << endmsg;
  //  theQaStream <<"no of fitted events with scaling factor: " << integralFitWeight*scaleFactor << "\n";

  std::map<std::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess >::iterator itMassMap;
  for(itMassMap= _massFitHistMap.begin(); itMassMap!= _massFitHistMap.end(); ++itMassMap){
    itMassMap->second->Scale(scaleFactor);
    if (itMassMap == _massFitHistMap.begin()) {
      //       InfoMsg << "No of fit events: " << itMassMap->second->Integral();
      //       theQaStream << "No of fit events: " << itMassMap->second->Integral() << "\n";
    }
  }
  //  theQaStream.close();

  std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess >::iterator itAngleMap;
  for(itAngleMap= _angleFitHistMap.begin(); itAngleMap!=_angleFitHistMap.end(); ++itAngleMap){
    std::vector<TH1F*>::iterator itTH1F;
    for(itTH1F=itAngleMap->second.begin(); itTH1F!=itAngleMap->second.end(); ++itTH1F){
      (*itTH1F)->Scale(scaleFactor);
    }
  }

  for(itAngleMap= _angleGJFitHistMap.begin(); itAngleMap!=_angleGJFitHistMap.end(); ++itAngleMap){
    std::vector<TH1F*>::iterator itTH1F;
    for(itTH1F=itAngleMap->second.begin(); itTH1F!=itAngleMap->second.end(); ++itTH1F){
      (*itTH1F)->Scale(scaleFactor);
    }
  }

  std::map<std::shared_ptr<angleHistData2D>, std::vector<TH2F*>, pawian::Collection::SharedPtrLess >::iterator itAngleMap2D;
  for(itAngleMap2D= _angleFitHistMap2D.begin(); itAngleMap2D!=_angleFitHistMap2D.end(); ++itAngleMap2D){
    std::vector<TH2F*>::iterator itTH2F;
    for(itTH2F=itAngleMap2D->second.begin(); itTH2F!=itAngleMap2D->second.end(); ++itTH2F){
      (*itTH2F)->Scale(scaleFactor);
    }
  }


}

void RootHist::fillEvt(EvtData* theData, double weight, std::string evtType){
  TTree* theTree=0;
  if(evtType=="data"){
    theTree=_dataFourvecs;
    fillMassHists(theData, weight, _massDataHistMap);
    fillAngleHists(theData, weight, _angleDataHistMap, "heli");
    fillAngleHists(theData, weight, _angleGJDataHistMap, "GottfriedJackson");
    fillAngleHists2D(theData, weight, _angleDataHistMap2D);
  }
  else if(evtType=="fit"){
    theTree=_fittedFourvecs;
    fillMassHists(theData, weight, _massFitHistMap);
    fillAngleHists(theData, weight, _angleFitHistMap, "heli");
    fillAngleHists(theData, weight, _angleGJFitHistMap, "GottfriedJackson");
    fillAngleHists2D(theData, weight, _angleFitHistMap2D);
  }
  else if(evtType=="mc"){
    fillMassHists(theData, weight, _massMcHistMap);
    fillAngleHists(theData, weight, _angleMcHistMap, "heli");
    fillAngleHists(theData, weight, _angleGJMcHistMap, "GottfriedJackson");
    fillAngleHists2D(theData, weight, _angleMcHistMap2D);
  }
  else if(evtType=="truthWoWeight"){
    fillMassHists(theData, weight, _massTruthHistMap);
    fillAngleHists(theData, weight, _angleTruthHistMap, "heli");
    fillAngleHists(theData, weight, _angleGJTruthHistMap, "GottfriedJackson");
    //    fillAngleHists2D(theData, weight, _angleMcHistMap2D);
  }
  else if(evtType=="truthWWeight"){
    theTree=_truthFourvecs;
    fillMassHists(theData, weight, _massTruthFitHistMap);
    fillAngleHists(theData, weight, _angleTruthFitHistMap, "heli");
    fillAngleHists(theData, weight, _angleGJTruthFitHistMap, "GottfriedJackson");
    //    fillAngleHists2D(theData, weight, _angleMcHistMap2D);
  }  
  if(evtType=="data" || evtType=="fit" || evtType=="truthWWeight"){
    for(auto fsIt = _fsParticles.begin(); fsIt != _fsParticles.end(); ++fsIt){
      std::string particleName = (*fsIt)->name();
      Vector4<double> tmp4vec=theData->FourVecsId.at(IdStringMapRegistry::instance()->stringId(particleName));
      _fourVecMap[particleName]->SetPxPyPzE(tmp4vec.X(), tmp4vec.Y(), tmp4vec.Z(), tmp4vec.E());
    }
    _weightToWrite = weight;
    theTree->Fill();
  }
}


void RootHist::fillMassHists(EvtData* theData, double weight, std::map<std::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess >& toFill){

  std::map<std::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess >::iterator it;
  for(it= toFill.begin(); it!= toFill.end(); ++it){
    //get relevant 4 vecs
    Vector4<double> combined4Vec(0.,0.,0.,0.);

    std::vector<std::string> fspNames=it->first->_fspNames;
    std::vector<std::string>::iterator itStr;
    for(itStr=fspNames.begin(); itStr!=fspNames.end(); ++itStr){
      Vector4<double> tmp4vec=theData->FourVecsId.at(IdStringMapRegistry::instance()->stringId(*itStr));
      combined4Vec+=tmp4vec;
    }
    it->second->Fill(combined4Vec.M(), weight);
  }

}

void RootHist::fillAngleHists(EvtData* theData, double weight, std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess >& toFill, std::string frame){

  std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess >::iterator it;
  for(it= toFill.begin(); it!= toFill.end(); ++it){
    short nBodyDecay = it->first->_nBodyDecay;

    Vector4<double> combinedDec4Vec(0.,0.,0.,0.);
    Vector4<double> combinedDec4Vec2(0.,0.,0.,0.);

    Vector4<double> combinedMother4Vec(0.,0.,0.,0.);
    std::string stringAll="all";
    Vector4<double> all4Vec=theData->FourVecsId.at(IdStringMapRegistry::instance()->stringId(stringAll));
    std::vector<std::string> decNames=it->first->_decPNames;
    std::vector<std::string> decNames2=it->first->_decPNames2;
    std::vector<std::string>::iterator itStr;
    for(itStr=decNames.begin(); itStr!=decNames.end(); ++itStr){
      Vector4<double> tmp4vec=theData->FourVecsId.at(IdStringMapRegistry::instance()->stringId(*itStr));
      combinedDec4Vec+=tmp4vec;
    }
    for(itStr=decNames2.begin(); itStr!=decNames2.end(); ++itStr){
      Vector4<double> tmp4vec=theData->FourVecsId.at(IdStringMapRegistry::instance()->stringId(*itStr));
      combinedDec4Vec2+=tmp4vec;
    }

    std::vector<std::string> motherNames=it->first->_motherPNames;

    for(itStr=motherNames.begin(); itStr!=motherNames.end(); ++itStr){
      Vector4<double> tmp4vec=theData->FourVecsId.at(IdStringMapRegistry::instance()->stringId(*itStr));
      combinedMother4Vec+=tmp4vec;
    }

    Vector4<double>  motherRef4Vec(0.,0.,0.,1.);
    Vector4<double>  result4Vec(0.,0.,0.,0.);
    Vector4<double>  result4Vec2(0.,0.,0.,0.);

    if (frame=="heli"){
      result4Vec=KinUtils::heliVec(motherRef4Vec, all4Vec, combinedMother4Vec, combinedDec4Vec);
      if(nBodyDecay==3) result4Vec2=KinUtils::heliVec(motherRef4Vec, all4Vec, combinedMother4Vec, combinedDec4Vec2);
    }
    else if (frame=="GottfriedJackson"){
      result4Vec=KinUtils::gottfriedJacksonVec(motherRef4Vec, all4Vec, combinedMother4Vec, combinedDec4Vec);
      if(nBodyDecay==3) result4Vec2=KinUtils::gottfriedJacksonVec(motherRef4Vec, all4Vec, combinedMother4Vec, combinedDec4Vec2);
    }
    else {
      Alert << "transformation into the frame " << frame << " not supported!!!" << endmsg;
      exit(0);
    }
    // if( fabs(all4Vec.E()-combinedMother4Vec.E()) < 1e-5
    // 	&& fabs(all4Vec.Px()-combinedMother4Vec.Px()) < 1e-5
    // 	&& fabs(all4Vec.Py()-combinedMother4Vec.Py()) < 1e-5
    // 	&& fabs(all4Vec.Pz()-combinedMother4Vec.Pz()) < 1e-5  ){
    // result4Vec=combinedDec4Vec;
    // result4Vec.Boost(all4Vec);
    // }
    // else{
    //   result4Vec=helicityVec(all4Vec, combinedMother4Vec, combinedDec4Vec);
    //    if(nBodyDecay==3)
    // 	  result4Vec2=helicityVec(all4Vec, combinedMother4Vec, combinedDec4Vec2);
    // }

    if(nBodyDecay == 2){
      it->second[0]->Fill( result4Vec.CosTheta(), weight);
      it->second[1]->Fill( result4Vec.Phi(), weight);
    }
    else if(nBodyDecay == 3){
      Vector4<double> result4VecPart1 = result4Vec;
      Vector4<double> result4VecPart2 = result4Vec2;
      //      std::cout << "\nresult4Vec: " << result4Vec << std::endl;
      //      std::cout << "result4Vec2: " << result4Vec2 << std::endl;
      Vector4<double> normVec(0,
      			     result4VecPart1.Y()*result4VecPart2.Z()-result4VecPart1.Z()*result4VecPart2.Y(),
      			     result4VecPart1.Z()*result4VecPart2.X()-result4VecPart1.X()*result4VecPart2.Z(),
      			     result4VecPart1.X()*result4VecPart2.Y()-result4VecPart1.Y()*result4VecPart2.X());

      it->second[0]->Fill( normVec.CosTheta(), weight);
      it->second[1]->Fill( normVec.Phi(), weight);
      //      if (frame=="heli"){
	Vector4<double> combinedMotherHeli4Vec(combinedMother4Vec.M(), 0., 0., 0.);
	Vector4<double> result4VecPart3(combinedMotherHeli4Vec.E()-result4VecPart1.E()-result4VecPart2.E()
					,combinedMotherHeli4Vec.Px()-result4VecPart1.Px()-result4VecPart2.Px()
					,combinedMotherHeli4Vec.Py()-result4VecPart1.Py()-result4VecPart2.Py()
					,combinedMotherHeli4Vec.Pz()-result4VecPart1.Pz()-result4VecPart2.Pz());
	double theQ=result4VecPart1.E()-result4VecPart1.M()+result4VecPart2.E()-result4VecPart2.M()+result4VecPart3.E()-result4VecPart3.M();
	double lambdaNorm=theQ*theQ*(theQ*theQ/108.+result4VecPart1.M()*theQ/9.+result4VecPart1.M()*result4VecPart1.M()/3.);
	double lambda=normVec.P()*normVec.P()/lambdaNorm;
	it->second[2]->Fill( lambda, weight);
	// }
    }
  }
}

void RootHist::fillAngleHists2D(EvtData* theData, double weight, std::map<std::shared_ptr<angleHistData2D>, std::vector<TH2F*>, pawian::Collection::SharedPtrLess >& toFill){

  std::map<std::shared_ptr<angleHistData2D>, std::vector<TH2F*>, pawian::Collection::SharedPtrLess >::iterator it;
  for(it= toFill.begin(); it!= toFill.end(); ++it){

    Vector4<double> combinedDec4Vec(0.,0.,0.,0.);
    Vector4<double> combinedMother4Vec(0.,0.,0.,0.);
    Vector4<double> combinedDec4Vec_2(0.,0.,0.,0.);
    Vector4<double> combinedMother4Vec_2(0.,0.,0.,0.);
    std::string stringAll="all";
    Vector4<double> all4Vec=theData->FourVecsId.at(IdStringMapRegistry::instance()->stringId(stringAll));

    std::vector<std::string> decNames=it->first->_decPNames;
    std::vector<std::string> decNames_2=it->first->_decPNames_2;
    std::vector<std::string> motherNames=it->first->_motherPNames;
    std::vector<std::string> motherNames_2=it->first->_motherPNames_2;

    std::vector<std::string>::iterator itStr;
    for(itStr=decNames.begin(); itStr!=decNames.end(); ++itStr){
      Vector4<double> tmp4vec=theData->FourVecsId.at(IdStringMapRegistry::instance()->stringId(*itStr));
      combinedDec4Vec+=tmp4vec;
    }

    for(itStr=motherNames.begin(); itStr!=motherNames.end(); ++itStr){
      Vector4<double> tmp4vec=theData->FourVecsId.at(IdStringMapRegistry::instance()->stringId(*itStr));
      combinedMother4Vec+=tmp4vec;
    }

    Vector4<double>  result4Vec(0.,0.,0.,0.);
    if( fabs(all4Vec.E()-combinedMother4Vec.E()) < 1e-5
	&& fabs(all4Vec.Px()-combinedMother4Vec.Px()) < 1e-5
	&& fabs(all4Vec.Py()-combinedMother4Vec.Py()) < 1e-5
	&& fabs(all4Vec.Pz()-combinedMother4Vec.Pz()) < 1e-5  ){
      result4Vec=combinedDec4Vec;
      result4Vec.Boost(all4Vec);
    }
    else{
       result4Vec=helicityVec(all4Vec, combinedMother4Vec, combinedDec4Vec);
    }

    for(itStr=decNames_2.begin(); itStr!=decNames_2.end(); ++itStr){
      Vector4<double> tmp4vec=theData->FourVecsId.at(IdStringMapRegistry::instance()->stringId(*itStr));
      combinedDec4Vec_2+=tmp4vec;
    }

    for(itStr=motherNames_2.begin(); itStr!=motherNames_2.end(); ++itStr){
      Vector4<double> tmp4vec=theData->FourVecsId.at(IdStringMapRegistry::instance()->stringId(*itStr));
      combinedMother4Vec_2+=tmp4vec;
    }

    Vector4<double>  result4Vec_2(0.,0.,0.,0.);
    if( fabs(all4Vec.E()-combinedMother4Vec_2.E()) < 1e-5
        && fabs(all4Vec.Px()-combinedMother4Vec_2.Px()) < 1e-5
        && fabs(all4Vec.Py()-combinedMother4Vec_2.Py()) < 1e-5
        && fabs(all4Vec.Pz()-combinedMother4Vec_2.Pz()) < 1e-5  ){
      result4Vec_2=combinedDec4Vec_2;
      result4Vec_2.Boost(all4Vec);
    }
    else{
      result4Vec_2=helicityVec(all4Vec, combinedMother4Vec_2, combinedDec4Vec_2);
    }

    it->second[0]->Fill( result4Vec.CosTheta(), result4Vec_2.CosTheta(),weight);
    it->second[1]->Fill( result4Vec.Phi(), result4Vec_2.Phi(),weight);
  }
}


void RootHist::scaleFitHists(double scaleFactor){

  std::map<std::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess >::iterator itMassMap;
  for(itMassMap= _massFitHistMap.begin(); itMassMap!= _massFitHistMap.end(); ++itMassMap){
    itMassMap->second->Scale(scaleFactor);
  }

  std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess >::iterator itAngleMap;
  for(itAngleMap= _angleFitHistMap.begin(); itAngleMap!=_angleFitHistMap.end(); ++itAngleMap){
    std::vector<TH1F*>::iterator itTH1F;
    for(itTH1F=itAngleMap->second.begin(); itTH1F!=itAngleMap->second.end(); ++itTH1F){
      (*itTH1F)->Scale(scaleFactor);
    }
  }

  for(itAngleMap= _angleGJFitHistMap.begin(); itAngleMap!=_angleGJFitHistMap.end(); ++itAngleMap){
    std::vector<TH1F*>::iterator itTH1F;
    for(itTH1F=itAngleMap->second.begin(); itTH1F!=itAngleMap->second.end(); ++itTH1F){
      (*itTH1F)->Scale(scaleFactor);
    }
  }
    
  std::map<std::shared_ptr<angleHistData2D>, std::vector<TH2F*>, pawian::Collection::SharedPtrLess >::iterator itAngleMap2D;
  for(itAngleMap2D= _angleFitHistMap2D.begin(); itAngleMap2D!=_angleFitHistMap2D.end(); ++itAngleMap2D){
    std::vector<TH2F*>::iterator itTH2F;
    for(itTH2F=itAngleMap2D->second.begin(); itTH2F!=itAngleMap2D->second.end(); ++itTH2F){
      (*itTH2F)->Scale(scaleFactor);
    }
  }
  
  if(_fillTruths){
    for(itMassMap=_massTruthFitHistMap.begin(); itMassMap!= _massTruthFitHistMap.end(); ++itMassMap){
      itMassMap->second->Scale(scaleFactor);
    }
    
    std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess >::iterator itAngleMap;
    for(itAngleMap= _angleTruthFitHistMap.begin(); itAngleMap!=_angleTruthFitHistMap.end(); ++itAngleMap){
      std::vector<TH1F*>::iterator itTH1F;
      for(itTH1F=itAngleMap->second.begin(); itTH1F!=itAngleMap->second.end(); ++itTH1F){
	(*itTH1F)->Scale(scaleFactor);
      }
    }
    
    for(itAngleMap= _angleGJTruthFitHistMap.begin(); itAngleMap!=_angleGJTruthFitHistMap.end(); ++itAngleMap){
      std::vector<TH1F*>::iterator itTH1F;
      for(itTH1F=itAngleMap->second.begin(); itTH1F!=itAngleMap->second.end(); ++itTH1F){
	(*itTH1F)->Scale(scaleFactor);
      }
    }

    // std::map<std::shared_ptr<angleHistData2D>, std::vector<TH2F*>, pawian::Collection::SharedPtrLess >::iterator itAngleMap2D;
    // for(itAngleMap2D= _angleFitHistMap2D.begin(); itAngleMap2D!=_angleFitHistMap2D.end(); ++itAngleMap2D){
    //   std::vector<TH2F*>::iterator itTH2F;
    //   for(itTH2F=itAngleMap2D->second.begin(); itTH2F!=itAngleMap2D->second.end(); ++itTH2F){
    // 	(*itTH2F)->Scale(scaleFactor);
    //   }
    // }

  }

}

void RootHist::initAngleHists(std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess >& theMap, std::shared_ptr<angleHistData> theHistData, std::string dataType, std::string systemType){
    std::string tmpBaseName= theHistData->_name;
    boost::replace_all(tmpBaseName,"+","p");
    boost::replace_all(tmpBaseName,"-","m");

    //helicity data
    std::string histThetaName=dataType+"Theta"+systemType+"_"+tmpBaseName;
    std::string histPhiName=dataType+"Phi"+systemType+"_"+tmpBaseName;
    std::string histThetaDescription = "cos#Theta " +theHistData->_name +systemType+dataType;
    std::string histPhiDescription = "#phi " +theHistData->_name +systemType+dataType;

    TH1F* currentThetaAngleDataHist=new TH1F(histThetaName.c_str(), histThetaDescription.c_str(), 100., -1., 1.);
    TH1F* currentPhiAngleDataHist=new TH1F(histPhiName.c_str(), histPhiDescription.c_str(), 100., -3.14159, 3.14159);
    currentThetaAngleDataHist->Sumw2();
    currentPhiAngleDataHist->Sumw2();
    theMap[theHistData].push_back(currentThetaAngleDataHist);
    theMap[theHistData].push_back(currentPhiAngleDataHist);

    if( theHistData->_nBodyDecay == 3){
      std::string histLambdaName="Lambda"+dataType+systemType+tmpBaseName;
      std::string histLambdaDescription = "#Lambda(" +theHistData->_name+dataType+systemType;
      TH1F* currentLambdaDataHist=new TH1F(histLambdaName.c_str(), histLambdaDescription.c_str(), 100., 0., 1.);
      currentLambdaDataHist->Sumw2();
      theMap[theHistData].push_back(currentLambdaDataHist);
    }
}


