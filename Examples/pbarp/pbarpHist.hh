#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <sstream>
#include <boost/shared_ptr.hpp>

#include "TROOT.h"


#include "PwaUtils/DataUtils.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/FitParamsBaseNew.hh"


//#include "TFile.h"
//#include "TH1F.h"
//#include "TH2F.h"

//class TFile;
//class TH2F;
//class TH1F;

class AbsLhNew;
class EvtDataNew;


struct massHistData {
  massHistData(std::vector<std::string>& fspNames) :
    _name("")
    ,_fspNames(fspNames)
  {
     std::vector<std::string>::iterator it;
    for(it=fspNames.begin(); it!=fspNames.end(); ++it){
      _name+=(*it);
    }
  }

  std::string _name;
  std::vector<std::string> _fspNames; 

  virtual bool operator==(const massHistData& compare) const {
    bool result=false;
    if ( _name==compare._name) result=true;
    return result;
  }

 virtual bool operator<(const massHistData& compare) const {
   bool result=false;
   if(_name < compare._name) result=true; 
    return result; 
  }  
};



class TFile;
//class TH2F;
class TH1F;
// class TNtuple;


class pbarpHist {

public:
  pbarpHist(boost::shared_ptr<AbsLhNew>, fitParamsNew& theFitParams);
  virtual ~pbarpHist();

protected:

private:
  void fillMassHists(EvtDataNew* theData, double weight);

 TFile* _theTFile;

 std::map<boost::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess > _massHistMap;

  void initRootStuff();
};

