#include <getopt.h>
#include <fstream>
#include <algorithm>

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
      fillMassHists((*it), 1.);
      ++it;
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
    std::string histName="data"+tmpMassHistData->_name;
    TH1F* currentMassHist=new TH1F(histName.c_str(), histName.c_str(), 100., 0., 2.);
    _massHistMap[tmpMassHistData]=currentMassHist;
  }

}

void pbarpHist::fillMassHists(EvtDataNew* theData, double weight){

  std::map<boost::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess >::iterator it;
  for(it= _massHistMap.begin(); it!= _massHistMap.end(); ++it){
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
