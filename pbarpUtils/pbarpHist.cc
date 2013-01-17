// pbarpHist class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>

#include "pbarpUtils/pbarpHist.hh"
#include "pbarpUtils/pbarpEnv.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/EvtDataBaseList.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
//#include "TMath.h"

pbarpHist::pbarpHist(boost::shared_ptr<AbsLh> theLh, fitParams& theFitParams) 
{
  if(0==theLh){
    Alert <<"AbsLh* is a 0 pointer !!!!" ;  // << endmsg;
    exit(1);
  }

  boost::shared_ptr<const EvtDataBaseList> theEvtList=theLh->getEventList();
  const std::vector<EvtData*> dataList=theEvtList->getDataVecs();

  initRootStuff();
  
  std::vector<EvtData*>::const_iterator it=dataList.begin();
  while(it!=dataList.end())
    {
      double weight = (*it)->evtWeight;
      fillMassHists((*it), weight, _massDataHistMap);
      fillAngleHists((*it), weight, _angleDataHistMap);
      ++it;
    }

  const std::vector<EvtData*> mcList=theEvtList->getMcVecs();  
  it=mcList.begin();
  while(it!=mcList.end())
    {
      double evtWeight = (*it)->evtWeight;
      fillMassHists((*it), evtWeight, _massMcHistMap);
      fillAngleHists((*it), evtWeight, _angleMcHistMap);

      double fitWeight= theLh->calcEvtIntensity( (*it), theFitParams );
      fillMassHists((*it), evtWeight*fitWeight, _massFitHistMap);
      fillAngleHists((*it), evtWeight*fitWeight, _angleFitHistMap);
      ++it;
    }

  double integralDataWoWeight=(double) theEvtList->getDataVecs().size();
  Info <<"No of data events without weight " << integralDataWoWeight ;  // << endmsg;

  double integralDataWWeight=(double) theEvtList->NoOfWeightedDataEvts();
  Info <<"No of data events with weight " << integralDataWWeight ;  // << endmsg;


  double integralMC=(double) theEvtList->NoOfWeightedMcEvts();
  Info <<"No of MC events " << integralMC ;  // << endmsg; 

  Info <<"scaling factor  " << integralDataWWeight/integralMC ;  // << endmsg;

  double scaleFactor = integralDataWWeight/integralMC;

  std::map<boost::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess >::iterator itMassMap;
  for(itMassMap= _massFitHistMap.begin(); itMassMap!= _massFitHistMap.end(); ++itMassMap){
    itMassMap->second->Scale(scaleFactor);
  }

  std::map<boost::shared_ptr<angleHistData>, std::pair<TH1F*, TH1F*>, pawian::Collection::SharedPtrLess >::iterator itAngleMap;
  for(itAngleMap= _angleFitHistMap.begin(); itAngleMap!=_angleFitHistMap.end(); ++itAngleMap){
    itAngleMap->second.first->Scale(scaleFactor);
    itAngleMap->second.second->Scale(scaleFactor);
  }
}

pbarpHist::~pbarpHist(){
  _theTFile->Write();
  _theTFile->Close();
}

void pbarpHist::initRootStuff(){
  std::ostringstream rootFileName;
  rootFileName << "./pawianHists" << pbarpEnv::instance()->outputFileNameSuffix() << ".root";
  _theTFile=new TFile(rootFileName.str().c_str(),"recreate");

  std::vector<std::vector<std::string> > histMassNameVec=pbarpEnv::instance()->histMassSystems();
  std::vector<std::vector<std::string> >::iterator itVecStr;
  for(itVecStr=histMassNameVec.begin(); itVecStr!=histMassNameVec.end(); ++itVecStr){
    boost::shared_ptr<massHistData> tmpMassHistData(new massHistData(*itVecStr));
    std::string tmpBaseName=tmpMassHistData->_name;
    boost::replace_all(tmpBaseName,"+","p");
    boost::replace_all(tmpBaseName,"-","m");
    std::string histName="data"+tmpBaseName;
    std::string histDescription = "M("+tmpMassHistData->_name+") (data)";

    double pMass = pbarpEnv::instance()->particleTable()->particle("proton")->mass();
    double pbarMom = pbarpEnv::instance()->pbarMomentum();
    double massMin = 0;
    double massMax = sqrt(pow(sqrt(pMass*pMass + pbarMom*pbarMom) + pMass, 2) - pbarMom*pbarMom);

    std::vector<std::string> fspNames=tmpMassHistData->_fspNames;
    std::vector<Particle*> allFsp =  pbarpEnv::instance()->finalStateParticles();
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
  }

  std::vector<boost::shared_ptr<angleHistData> > angleHistDataVec=pbarpEnv::instance()->angleHistDataVec();

  std::vector<boost::shared_ptr<angleHistData> >::iterator itAngleVec;
  for (itAngleVec=angleHistDataVec.begin(); itAngleVec!=angleHistDataVec.end(); ++itAngleVec){
    std::string tmpBaseName= (*itAngleVec)->_name;
    boost::replace_all(tmpBaseName,"+","p");
    boost::replace_all(tmpBaseName,"-","m");
    std::string histThetaName="dataTheta"+tmpBaseName;
    std::string histPhiName="dataPhi"+tmpBaseName;
    std::string histThetaDescription = "cos#Theta(" +(*itAngleVec)->_name + ") (data)";
    std::string histPhiDescription = "#phi(" +(*itAngleVec)->_name + ") (data)";

    TH1F* currentThetaAngleDataHist=new TH1F(histThetaName.c_str(), histThetaDescription.c_str(), 100., -1., 1.);
    TH1F* currentPhiAngleDataHist=new TH1F(histPhiName.c_str(), histPhiDescription.c_str(), 100., -3.14159, 3.14159);
    currentThetaAngleDataHist->Sumw2();
    currentPhiAngleDataHist->Sumw2();
    _angleDataHistMap[*itAngleVec]=std::pair<TH1F*, TH1F*>(currentThetaAngleDataHist, currentPhiAngleDataHist);

    histThetaName="MCTheta"+tmpBaseName;
    histPhiName="MCPhi"+tmpBaseName;
    histThetaDescription = "cos#Theta(" +(*itAngleVec)->_name + ") (MC)";
    histPhiDescription = "#phi(" +(*itAngleVec)->_name + ") (MC)";

    TH1F* currentThetaAngleMcHist=new TH1F(histThetaName.c_str(), histThetaDescription.c_str(), 100., -1., 1.);
    TH1F* currentPhiAngleMcHist=new TH1F(histPhiName.c_str(), histPhiDescription.c_str(), 100., -3.14159, 3.14159);
    currentThetaAngleMcHist->Sumw2();
    currentPhiAngleMcHist->Sumw2();
    _angleMcHistMap[*itAngleVec]=std::pair<TH1F*, TH1F*>(currentThetaAngleMcHist, currentPhiAngleMcHist);

    histThetaName="FitTheta"+tmpBaseName;
    histPhiName="FitPhi"+tmpBaseName;
    histThetaDescription = "cos#Theta(" +(*itAngleVec)->_name + ") (fit)";
    histPhiDescription = "#phi(" +(*itAngleVec)->_name + ") (fit)";

    TH1F* currentThetaAngleFitHist=new TH1F(histThetaName.c_str(), histThetaDescription.c_str(), 100., -1., 1.);
    TH1F* currentPhiAngleFitHist=new TH1F(histPhiName.c_str(), histPhiDescription.c_str(), 100., -3.14159, 3.14159);
    currentThetaAngleFitHist->Sumw2();
    currentPhiAngleFitHist->Sumw2();
    _angleFitHistMap[*itAngleVec]=std::pair<TH1F*, TH1F*>(currentThetaAngleFitHist, currentPhiAngleFitHist);
  } 

 std::map<boost::shared_ptr<angleHistData>, TH1F*, pawian::Collection::SharedPtrLess > _angleDataHistMap;
}

void pbarpHist::fillMassHists(EvtData* theData, double weight, std::map<boost::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess >& toFill){

  std::map<boost::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess >::iterator it;
  for(it= toFill.begin(); it!= toFill.end(); ++it){
    //get relevant 4 vecs
    Vector4<double> combined4Vec(0.,0.,0.,0.);

    std::vector<std::string> fspNames=it->first->_fspNames;
    std::vector<std::string>::iterator itStr;
    for(itStr=fspNames.begin(); itStr!=fspNames.end(); ++itStr){
      Vector4<double> tmp4vec=theData->FourVecsString[*itStr];
      combined4Vec+=tmp4vec;
    }
    it->second->Fill(combined4Vec.M(), weight); 
  }
  
}

void pbarpHist::fillAngleHists(EvtData* theData, double weight, std::map<boost::shared_ptr<angleHistData>, std::pair<TH1F*, TH1F*>, pawian::Collection::SharedPtrLess >& toFill){

  std::map<boost::shared_ptr<angleHistData>, std::pair<TH1F*, TH1F*>, pawian::Collection::SharedPtrLess >::iterator it;
  for(it= toFill.begin(); it!= toFill.end(); ++it){
    Vector4<double> combinedDec4Vec(0.,0.,0.,0.);
    Vector4<double> combinedMother4Vec(0.,0.,0.,0.);
    Vector4<double> all4Vec=theData->FourVecsString["all"];

    std::vector<std::string> decNames=it->first->_decPNames;
    std::vector<std::string>::iterator itStr;
    for(itStr=decNames.begin(); itStr!=decNames.end(); ++itStr){
      Vector4<double> tmp4vec=theData->FourVecsString[*itStr];
      combinedDec4Vec+=tmp4vec;
    }

    std::vector<std::string> motherNames=it->first->_motherPNames;

    for(itStr=motherNames.begin(); itStr!=motherNames.end(); ++itStr){
      Vector4<double> tmp4vec=theData->FourVecsString[*itStr];
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
    else result4Vec=helicityVec(all4Vec, combinedMother4Vec, combinedDec4Vec);

    it->second.first->Fill( result4Vec.CosTheta(), weight);  
    it->second.second->Fill( result4Vec.Phi(), weight);
  }
}
