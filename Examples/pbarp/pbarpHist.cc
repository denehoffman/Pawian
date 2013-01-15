#include <getopt.h>
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>

#include "Examples/pbarp/pbarpHist.hh"
#include "Examples/pbarp/pbarpEnv.hh"
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
      fillThetaHists((*it), weight, _thetaDataHistMap);
      ++it;
    }

  const std::vector<EvtData*> mcList=theEvtList->getMcVecs();  
  it=mcList.begin();
  while(it!=mcList.end())
    {
      double evtWeight = (*it)->evtWeight;
      fillMassHists((*it), evtWeight, _massMcHistMap);
      fillThetaHists((*it), evtWeight, _thetaMcHistMap);

      double fitWeight= theLh->calcEvtIntensity( (*it), theFitParams );
      fillMassHists((*it), evtWeight*fitWeight, _massFitHistMap);
      fillThetaHists((*it), evtWeight*fitWeight, _thetaFitHistMap);
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

  std::map<boost::shared_ptr<angleHistData>, TH1F*, pawian::Collection::SharedPtrLess >::iterator itAngleMap;
  for(itAngleMap= _thetaFitHistMap.begin(); itAngleMap!=_thetaFitHistMap.end(); ++itAngleMap){
    itAngleMap->second->Scale(scaleFactor);
  }
}

pbarpHist::~pbarpHist(){
  _theTFile->Write();
  _theTFile->Close();
}

void pbarpHist::initRootStuff(){
  std::string rootFileName="./pawianHists.root";
  _theTFile=new TFile(rootFileName.c_str(),"recreate");

  std::vector<std::vector<std::string> > histMassNameVec=pbarpEnv::instance()->histMassSystems();
  std::vector<std::vector<std::string> >::iterator itVecStr;
  for(itVecStr=histMassNameVec.begin(); itVecStr!=histMassNameVec.end(); ++itVecStr){
    boost::shared_ptr<massHistData> tmpMassHistData(new massHistData(*itVecStr));
    std::string tmpBaseName=tmpMassHistData->_name;
    boost::replace_all(tmpBaseName,"+","p");
    boost::replace_all(tmpBaseName,"-","m");
    std::string histName="data"+tmpBaseName;
    std::string histDescription = "M("+tmpMassHistData->_name+") (data)";

    double massMin=-0.05;
    double massMax=2.;
    std::vector<std::string> fspNames=tmpMassHistData->_fspNames;
    std::vector<std::string>::iterator itStr2;
    for(itStr2=fspNames.begin(); itStr2!=fspNames.end(); ++itStr2){
      Particle* currentParticle=pbarpEnv::instance()->particleTable()->particle(*itStr2);
      massMin+= currentParticle->mass(); 
      massMax+= currentParticle->mass(); 
    }

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
    std::string histName="data"+tmpBaseName;
    std::string histDescription = "cos#Theta(" +(*itAngleVec)->_name + ") (data)";
    TH1F* currentThetaAngleDataHist=new TH1F(histName.c_str(), histDescription.c_str(), 100., -1., 1.);
    currentThetaAngleDataHist->Sumw2();
    _thetaDataHistMap[*itAngleVec]=currentThetaAngleDataHist;

    histName="MC"+tmpBaseName;
    histDescription = "cos#Theta(" +(*itAngleVec)->_name + ") (MC)";
    TH1F* currentThetaAngleMcHist=new TH1F(histName.c_str(), histDescription.c_str(), 100., -1., 1.);
    currentThetaAngleMcHist->Sumw2();
    _thetaMcHistMap[*itAngleVec]=currentThetaAngleMcHist;

    histName="Fit"+tmpBaseName;
    histDescription = "cos#Theta(" +(*itAngleVec)->_name + ") (fit)";
    TH1F* currentThetaAngleFitHist=new TH1F(histName.c_str(), histDescription.c_str(), 100., -1., 1.);
    currentThetaAngleFitHist->Sumw2();
    _thetaFitHistMap[*itAngleVec]=currentThetaAngleFitHist;
  } 

 std::map<boost::shared_ptr<angleHistData>, TH1F*, pawian::Collection::SharedPtrLess > _thetaDataHistMap;
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

void pbarpHist::fillThetaHists(EvtData* theData, double weight, std::map<boost::shared_ptr<angleHistData>, TH1F*, pawian::Collection::SharedPtrLess >& toFill){

  std::map<boost::shared_ptr<angleHistData>, TH1F*, pawian::Collection::SharedPtrLess >::iterator it;
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

    it->second->Fill( result4Vec.CosTheta(), weight);     
  }
}
