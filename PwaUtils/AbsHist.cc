// AbsHist class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>

#include "PwaUtils/AbsHist.hh"
#include "PwaUtils/AbsEnv.hh"
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

AbsHist::AbsHist(AbsEnv* theEnv) :
  _absEnv(theEnv) 
{
  std::ostringstream rootFileName;
  rootFileName << "./pawianHists" << _absEnv->outputFileNameSuffix() << ".root";
  _theTFile=new TFile(rootFileName.str().c_str(),"recreate");

  std::vector<boost::shared_ptr<angleHistData> > angleHistDataVec=_absEnv->angleHistDataVec();

  std::vector<boost::shared_ptr<angleHistData> >::iterator itAngleVec;
  for (itAngleVec=angleHistDataVec.begin(); itAngleVec!=angleHistDataVec.end(); ++itAngleVec){
    std::string tmpBaseName= (*itAngleVec)->_name;
    boost::replace_all(tmpBaseName,"+","p");
    boost::replace_all(tmpBaseName,"-","m");
    std::string histThetaName="DataTheta"+tmpBaseName;
    std::string histPhiName="DataPhi"+tmpBaseName;
    std::string histThetaDescription = "cos#Theta(" +(*itAngleVec)->_name + ") (data)";
    std::string histPhiDescription = "#phi(" +(*itAngleVec)->_name + ") (data)";

    TH1F* currentThetaAngleDataHist=new TH1F(histThetaName.c_str(), histThetaDescription.c_str(), 100., -1., 1.);
    TH1F* currentPhiAngleDataHist=new TH1F(histPhiName.c_str(), histPhiDescription.c_str(), 100., -3.14159, 3.14159);
    currentThetaAngleDataHist->Sumw2();
    currentPhiAngleDataHist->Sumw2();
    _angleDataHistMap[*itAngleVec].push_back(currentThetaAngleDataHist);
    _angleDataHistMap[*itAngleVec].push_back(currentPhiAngleDataHist);

    histThetaName="MCTheta"+tmpBaseName;
    histPhiName="MCPhi"+tmpBaseName;
    histThetaDescription = "cos#Theta(" +(*itAngleVec)->_name + ") (MC)";
    histPhiDescription = "#phi(" +(*itAngleVec)->_name + ") (MC)";

    TH1F* currentThetaAngleMcHist=new TH1F(histThetaName.c_str(), histThetaDescription.c_str(), 100., -1., 1.);
    TH1F* currentPhiAngleMcHist=new TH1F(histPhiName.c_str(), histPhiDescription.c_str(), 100., -3.14159, 3.14159);
    currentThetaAngleMcHist->Sumw2();
    currentPhiAngleMcHist->Sumw2();
    _angleMcHistMap[*itAngleVec].push_back(currentThetaAngleMcHist);
    _angleMcHistMap[*itAngleVec].push_back(currentPhiAngleMcHist);

    histThetaName="FitTheta"+tmpBaseName;
    histPhiName="FitPhi"+tmpBaseName;
    histThetaDescription = "cos#Theta(" +(*itAngleVec)->_name + ") (fit)";
    histPhiDescription = "#phi(" +(*itAngleVec)->_name + ") (fit)";

    TH1F* currentThetaAngleFitHist=new TH1F(histThetaName.c_str(), histThetaDescription.c_str(), 100., -1., 1.);
    TH1F* currentPhiAngleFitHist=new TH1F(histPhiName.c_str(), histPhiDescription.c_str(), 100., -3.14159, 3.14159);
    currentThetaAngleFitHist->Sumw2();
    currentPhiAngleFitHist->Sumw2();
    _angleFitHistMap[*itAngleVec].push_back(currentThetaAngleFitHist);
    _angleFitHistMap[*itAngleVec].push_back(currentPhiAngleFitHist);

    if( (*itAngleVec)->_nBodyDecay == 3){
      std::string histLambdaName="DataLambda"+tmpBaseName;
      std::string histLambdaDescription = "#Lambda(" +(*itAngleVec)->_name + ") (data)";
      TH1F* currentLambdaDataHist=new TH1F(histLambdaName.c_str(), histLambdaDescription.c_str(), 100., 0., 1.);
      _angleDataHistMap[*itAngleVec].push_back(currentLambdaDataHist);

      histLambdaName="MCLambda"+tmpBaseName;
      histLambdaDescription = "#Lambda(" +(*itAngleVec)->_name + ") (MC)";
      TH1F* currentLambdaMcHist=new TH1F(histLambdaName.c_str(), histLambdaDescription.c_str(), 100., 0., 1.); 
      _angleMcHistMap[*itAngleVec].push_back(currentLambdaMcHist);
 
      histLambdaName="FitLambda"+tmpBaseName;
      histLambdaDescription = "#Lambda(" +(*itAngleVec)->_name + ") (fit)";
      TH1F* currentLambdaFitHist=new TH1F(histLambdaName.c_str(), histLambdaDescription.c_str(), 100., 0., 1.); 
      _angleFitHistMap[*itAngleVec].push_back(currentLambdaFitHist);
    }
  } 
  // 2D-Histograms

 std::vector<boost::shared_ptr<angleHistData2D> > angleHistDataVec2D=_absEnv->angleHistDataVec2D();

  std::vector<boost::shared_ptr<angleHistData2D> >::iterator itAngleVec2D;
  for (itAngleVec2D=angleHistDataVec2D.begin(); itAngleVec2D!=angleHistDataVec2D.end(); ++itAngleVec2D){
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
}

AbsHist::~AbsHist(){
  _theTFile->Write();
  _theTFile->Close();
}

void AbsHist::fillIt(boost::shared_ptr<AbsLh> theLh, fitParams& theFitParams){

  if(0==theLh){
    Alert <<"AbsLh* is a 0 pointer !!!!" ;  // << endmsg;
    exit(1);
  }

  boost::shared_ptr<const EvtDataBaseList> theEvtList=theLh->getEventList();
  const std::vector<EvtData*> dataList=theEvtList->getDataVecs();

  std::vector<EvtData*>::const_iterator it=dataList.begin();
  while(it!=dataList.end())
    {
      double weight = (*it)->evtWeight;
      fillMassHists((*it), weight, _massDataHistMap);
      fillAngleHists((*it), weight, _angleDataHistMap);
      fillAngleHists2D((*it), weight, _angleDataHistMap2D);
      ++it;
    }

  const std::vector<EvtData*> mcList=theEvtList->getMcVecs();  
  it=mcList.begin();
  while(it!=mcList.end())
    {
      double evtWeight = (*it)->evtWeight;
      fillMassHists((*it), evtWeight, _massMcHistMap);
      fillAngleHists((*it), evtWeight, _angleMcHistMap);
      fillAngleHists2D((*it), evtWeight, _angleMcHistMap2D);

      double fitWeight= theLh->calcEvtIntensity( (*it), theFitParams );
      fillMassHists((*it), evtWeight*fitWeight, _massFitHistMap);
      fillAngleHists((*it), evtWeight*fitWeight, _angleFitHistMap);
      fillAngleHists2D((*it), evtWeight*fitWeight, _angleFitHistMap2D);
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

  std::map<boost::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess >::iterator itAngleMap;
  for(itAngleMap= _angleFitHistMap.begin(); itAngleMap!=_angleFitHistMap.end(); ++itAngleMap){
    std::vector<TH1F*>::iterator itTH1F;
    for(itTH1F=itAngleMap->second.begin(); itTH1F!=itAngleMap->second.end(); ++itTH1F){  
      (*itTH1F)->Scale(scaleFactor);
    }
  }
  
  std::map<boost::shared_ptr<angleHistData2D>, std::vector<TH2F*>, pawian::Collection::SharedPtrLess >::iterator itAngleMap2D;
  for(itAngleMap2D= _angleFitHistMap2D.begin(); itAngleMap2D!=_angleFitHistMap2D.end(); ++itAngleMap2D){
    std::vector<TH2F*>::iterator itTH2F;
    for(itTH2F=itAngleMap2D->second.begin(); itTH2F!=itAngleMap2D->second.end(); ++itTH2F){  
      (*itTH2F)->Scale(scaleFactor);
    }
  }


}



void AbsHist::fillMassHists(EvtData* theData, double weight, std::map<boost::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess >& toFill){

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

void AbsHist::fillAngleHists(EvtData* theData, double weight, std::map<boost::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess >& toFill){

  std::map<boost::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess >::iterator it;
  for(it= toFill.begin(); it!= toFill.end(); ++it){
    short nBodyDecay = it->first->_nBodyDecay;

    Vector4<double> combinedDec4Vec(0.,0.,0.,0.);
    Vector4<double> combinedDec4Vec2(0.,0.,0.,0.);

    Vector4<double> combinedMother4Vec(0.,0.,0.,0.);
    Vector4<double> all4Vec=theData->FourVecsString["all"];

    std::vector<std::string> decNames=it->first->_decPNames;
    std::vector<std::string> decNames2=it->first->_decPNames2;
    std::vector<std::string>::iterator itStr;
    for(itStr=decNames.begin(); itStr!=decNames.end(); ++itStr){
      Vector4<double> tmp4vec=theData->FourVecsString[*itStr];
      combinedDec4Vec+=tmp4vec;
    }
    for(itStr=decNames2.begin(); itStr!=decNames2.end(); ++itStr){
      Vector4<double> tmp4vec=theData->FourVecsString[*itStr];
      combinedDec4Vec2+=tmp4vec;
    }

    std::vector<std::string> motherNames=it->first->_motherPNames;

    for(itStr=motherNames.begin(); itStr!=motherNames.end(); ++itStr){
      Vector4<double> tmp4vec=theData->FourVecsString[*itStr];
      combinedMother4Vec+=tmp4vec;
    }

    Vector4<double>  result4Vec(0.,0.,0.,0.);
    Vector4<double>  result4Vec2(0.,0.,0.,0.);

    if( fabs(all4Vec.E()-combinedMother4Vec.E()) < 1e-5
	&& fabs(all4Vec.Px()-combinedMother4Vec.Px()) < 1e-5 
	&& fabs(all4Vec.Py()-combinedMother4Vec.Py()) < 1e-5
	&& fabs(all4Vec.Pz()-combinedMother4Vec.Pz()) < 1e-5  ){
    result4Vec=combinedDec4Vec;
    result4Vec.Boost(all4Vec);
    } 
    else{
       result4Vec=helicityVec(all4Vec, combinedMother4Vec, combinedDec4Vec);
       if(nBodyDecay==3)
	  result4Vec2=helicityVec(all4Vec, combinedMother4Vec, combinedDec4Vec2);
    }

    if(nBodyDecay == 2){
      it->second[0]->Fill( result4Vec.CosTheta(), weight);
      it->second[1]->Fill( result4Vec.Phi(), weight);
    }
    else if(nBodyDecay == 3){
      Vector4<double> result4VecPart1 = result4Vec;
      Vector4<double> result4VecPart2 = result4Vec2;
      Vector4<float> normVec(0,
			     result4VecPart1.Y()*result4VecPart2.Z()-result4VecPart1.Z()*result4VecPart2.Y(),
			     result4VecPart1.Z()*result4VecPart2.X()-result4VecPart1.X()*result4VecPart2.Z(),
			     result4VecPart1.X()*result4VecPart2.Y()-result4VecPart1.Y()*result4VecPart2.X());

      it->second[0]->Fill( normVec.CosTheta(), weight);
      it->second[1]->Fill( normVec.Phi(), weight);

      Vector4<double> result4VecPart3(combinedMother4Vec.M()-result4VecPart1.E()-result4VecPart2.E()
				      ,-result4VecPart1.Px()-result4VecPart2.Px()
				      ,-result4VecPart1.Py()-result4VecPart2.Py()
				      ,-result4VecPart1.Pz()-result4VecPart2.Pz());

      double theQ=result4VecPart1.E()-result4VecPart1.M()+result4VecPart2.E()-result4VecPart2.M()+result4VecPart3.E()-result4VecPart3.M();
      double lambdaNorm=theQ*theQ*(theQ*theQ/108.+result4VecPart1.M()*theQ/9.+result4VecPart1.M()*result4VecPart1.M()/3.);
      double lambda=normVec.P()*normVec.P()/lambdaNorm;
      it->second[2]->Fill( lambda, weight);
    }
  }
}

void AbsHist::fillAngleHists2D(EvtData* theData, double weight, std::map<boost::shared_ptr<angleHistData2D>, std::vector<TH2F*>, pawian::Collection::SharedPtrLess >& toFill){

  std::map<boost::shared_ptr<angleHistData2D>, std::vector<TH2F*>, pawian::Collection::SharedPtrLess >::iterator it;
  for(it= toFill.begin(); it!= toFill.end(); ++it){

    Vector4<double> combinedDec4Vec(0.,0.,0.,0.);
    Vector4<double> combinedMother4Vec(0.,0.,0.,0.);
    Vector4<double> combinedDec4Vec_2(0.,0.,0.,0.);
    Vector4<double> combinedMother4Vec_2(0.,0.,0.,0.);
    Vector4<double> all4Vec=theData->FourVecsString["all"];

    std::vector<std::string> decNames=it->first->_decPNames;
    std::vector<std::string> decNames_2=it->first->_decPNames_2;
    std::vector<std::string> motherNames=it->first->_motherPNames;
    std::vector<std::string> motherNames_2=it->first->_motherPNames_2;

    std::vector<std::string>::iterator itStr;
    for(itStr=decNames.begin(); itStr!=decNames.end(); ++itStr){
      Vector4<double> tmp4vec=theData->FourVecsString[*itStr];
      combinedDec4Vec+=tmp4vec;
    }
    
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
    else{
       result4Vec=helicityVec(all4Vec, combinedMother4Vec, combinedDec4Vec);
    }
    
    for(itStr=decNames_2.begin(); itStr!=decNames_2.end(); ++itStr){
      Vector4<double> tmp4vec=theData->FourVecsString[*itStr];
      combinedDec4Vec_2+=tmp4vec;
    }

    for(itStr=motherNames_2.begin(); itStr!=motherNames_2.end(); ++itStr){
      Vector4<double> tmp4vec=theData->FourVecsString[*itStr];
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
