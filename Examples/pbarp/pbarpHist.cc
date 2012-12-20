#include <getopt.h>
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>

#include "Examples/pbarp/pbarpHist.hh"
#include "Examples/pbarp/pbarpEnv.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/AbsLhNew.hh"
#include "PwaUtils/EvtDataBaseListNew.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
//#include "TMath.h"

pbarpHist::pbarpHist(boost::shared_ptr<AbsLhNew> theLh, fitParamsNew& theFitParams) 
{
  if(0==theLh){
    Alert <<"AbsLh* is a 0 pointer !!!!" ;  // << endmsg;
    exit(1);
  }

  boost::shared_ptr<const EvtDataBaseListNew> theEvtList=theLh->getEventList();
  const std::vector<EvtDataNew*> dataList=theEvtList->getDataVecs();

  initRootStuff();
  
  std::vector<EvtDataNew*>::const_iterator it=dataList.begin();
  while(it!=dataList.end())
    {
      double weight = (*it)->evtWeight;
      fillMassHists((*it), weight, _massDataHistMap);
      ++it;
    }

  const std::vector<EvtDataNew*> mcList=theEvtList->getMcVecs();  
  it=mcList.begin();
  while(it!=mcList.end())
    {
      double evtWeight = (*it)->evtWeight;
      fillMassHists((*it), evtWeight, _massMcHistMap);

      double fitWeight= theLh->calcEvtIntensity( (*it), theFitParams );
      fillMassHists((*it), evtWeight*fitWeight, _massFitHistMap);
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
    TH1F* currentMassDataHist=new TH1F(histName.c_str(), histDescription.c_str(), 100., 0., 2.);
    currentMassDataHist->Sumw2();
    _massDataHistMap[tmpMassHistData]=currentMassDataHist;

    histName="MC"+tmpBaseName;
    histDescription = "M("+tmpMassHistData->_name+") (MC)";
    TH1F* currentMassMcHist=new TH1F(histName.c_str(), histDescription.c_str(), 100., 0., 2.);
    currentMassMcHist->Sumw2();
    _massMcHistMap[tmpMassHistData]=currentMassMcHist;

    histName="Fit"+tmpBaseName;
    histDescription = "M("+tmpMassHistData->_name+") (fit)";
    TH1F* currentMassFitHist=new TH1F(histName.c_str(), histDescription.c_str(), 100., 0., 2.);
    currentMassFitHist->Sumw2();
    _massFitHistMap[tmpMassHistData]=currentMassFitHist;
  }

}

void pbarpHist::fillMassHists(EvtDataNew* theData, double weight, std::map<boost::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess >& toFill){

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
