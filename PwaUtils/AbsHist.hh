// AbsHist class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

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
#include "PwaUtils/FitParamsBase.hh"


class AbsLh;
class EvtData;


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

struct angleHistData {
   angleHistData(std::vector<std::string>& motherPNames, std::vector<std::string>& decPNames, std::vector<std::string>& decPNames2, short nBodyDecay) :
     _nBodyDecay(nBodyDecay)
    ,_name("")
    ,_motherPNames(motherPNames)
    ,_decPNames(decPNames)
    ,_decPNames2(decPNames2)
  {

    std::vector<std::string>::iterator it;
    for(it=decPNames.begin(); it!=decPNames.end(); ++it){
      _name+=(*it);
    }

    if(nBodyDecay == 3)
       _name+="AND";

    for(it=decPNames2.begin(); it!=decPNames2.end(); ++it){
      _name+=(*it);
    }

    _name+="_Heli";

    for(it=motherPNames.begin(); it!=motherPNames.end(); ++it){
      _name+=(*it);
    }
  }

  short _nBodyDecay;
  std::string _name;
  std::vector<std::string> _motherPNames;
  std::vector<std::string> _decPNames;
  std::vector<std::string> _decPNames2;

  virtual bool operator==(const angleHistData& compare) const {
    bool result=false;
    if ( _name==compare._name) result=true;
    return result;
  }

 virtual bool operator<(const angleHistData& compare) const {
   bool result=false;
   if(_name < compare._name) result=true; 
    return result; 
  }  
};

class TFile;
//class TH2F;
class TH1F;


class AbsHist {

public:
  AbsHist();
  virtual ~AbsHist();

protected:

  void fillIt(boost::shared_ptr<AbsLh> theLh, fitParams& theFitParams);

  void fillMassHists(EvtData* theData, double weight, std::map<boost::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess >& toFill);

  void fillAngleHists(EvtData* theData, double weight, std::map<boost::shared_ptr<angleHistData>, std::pair<TH1F*, TH1F*>, pawian::Collection::SharedPtrLess >& toFill);

 TFile* _theTFile;

 std::map<boost::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess > _massDataHistMap;
 std::map<boost::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess > _massMcHistMap;
 std::map<boost::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess > _massFitHistMap;

 std::map<boost::shared_ptr<angleHistData>, std::pair<TH1F*, TH1F*>, pawian::Collection::SharedPtrLess > _angleDataHistMap;
 std::map<boost::shared_ptr<angleHistData>, std::pair<TH1F*, TH1F*>, pawian::Collection::SharedPtrLess > _angleMcHistMap;
 std::map<boost::shared_ptr<angleHistData>, std::pair<TH1F*, TH1F*>, pawian::Collection::SharedPtrLess > _angleFitHistMap;

 virtual void initRootStuff()=0;

private:

};

