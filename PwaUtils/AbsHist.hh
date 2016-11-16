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

// AbsHist class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <sstream>
#include <memory>

#include "TROOT.h"

#include "PwaUtils/DataUtils.hh"
#include "Utils/PawianCollectionUtils.hh"
//#include "FitParams/FitParColBase.hh"


class AbsLh;
class EvtData;
class TLorentzVector; 

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

    _name+="_From";

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

struct angleHistData2D {
  angleHistData2D(std::vector<std::string>& motherPNames, std::vector<std::string>& decPNames, std::vector<std::string>& motherPNames_2, std::vector<std::string>& decPNames_2 ) :
    _name("")
    ,_nameXAxis("")
    ,_nameYAxis("")
    ,_motherPNames(motherPNames)
    ,_decPNames(decPNames)
    ,_motherPNames_2(motherPNames_2)
    ,_decPNames_2(decPNames_2)
  {

    std::vector<std::string>::iterator it;
    for(it=decPNames.begin(); it!=decPNames.end(); ++it){
      _name+=(*it);
      _nameXAxis+=(*it);
    }

    _name+="_Heli";
    _nameXAxis+="_Heli";

    for(it=motherPNames.begin(); it!=motherPNames.end(); ++it){
      _name+=(*it);
      _nameXAxis+=(*it);
    }

    _name+="_VS_";

    for(it=decPNames_2.begin(); it!=decPNames_2.end(); ++it){
      _name+=(*it);
      _nameYAxis+=(*it);
    }

    _name+="_Heli";
    _nameYAxis+="_Heli";

    for(it=motherPNames_2.begin(); it!=motherPNames_2.end(); ++it){
      _name+=(*it);
      _nameYAxis+=(*it);
    }

  }

  std::string _name;
  std::string _nameXAxis;
  std::string _nameYAxis;
  std::vector<std::string> _motherPNames;
  std::vector<std::string> _decPNames;
  std::vector<std::string> _motherPNames_2;
  std::vector<std::string> _decPNames_2;

  virtual bool operator==(const angleHistData2D& compare) const {
    bool result=false;
    if ( _name==compare._name) result=true;
    return result;
  }

 virtual bool operator<(const angleHistData2D& compare) const {
   bool result=false;
   if(_name < compare._name) result=true;
    return result;
  }
};


class TFile;
class TH2F;
class TH1F;
class TTree;
class Particle;
class AbsPawianParameters;

class AbsHist {

public:
  AbsHist(std::string additionalSuffix = "");
  virtual ~AbsHist();
  void fillEvt(EvtData* theData, double weight, std::string evtType);
  void fillFromLhData(std::shared_ptr<AbsLh> theLh, std::shared_ptr<AbsPawianParameters> fitParams);
  void scaleFitHists(double scaleFactor);

protected:

  void fillMassHists(EvtData* theData, double weight, std::map<std::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess >& toFill);
  void fillAngleHists(EvtData* theData, double weight, std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess >& toFill, std::string frame="heli");
  void fillAngleHists2D(EvtData* theData, double weight, std::map<std::shared_ptr<angleHistData2D>, std::vector<TH2F*>, pawian::Collection::SharedPtrLess >& toFill);

 TFile* _theTFile;
 TTree* _dataFourvecs;
 TTree* _fittedFourvecs;

 std::map<std::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess > _massDataHistMap;
 std::map<std::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess > _massMcHistMap;
 std::map<std::shared_ptr<massHistData>, TH1F*, pawian::Collection::SharedPtrLess > _massFitHistMap;

 std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess > _angleDataHistMap;
 std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess > _angleMcHistMap;
 std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess > _angleFitHistMap;

 std::map<std::shared_ptr<angleHistData2D>, std::vector<TH2F*>, pawian::Collection::SharedPtrLess > _angleDataHistMap2D;
 std::map<std::shared_ptr<angleHistData2D>, std::vector<TH2F*>, pawian::Collection::SharedPtrLess > _angleMcHistMap2D;
 std::map<std::shared_ptr<angleHistData2D>, std::vector<TH2F*>, pawian::Collection::SharedPtrLess > _angleFitHistMap2D;

  //angles in the Gottfried Jackson system
  std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess > _angleGJDataHistMap;
 std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess > _angleGJMcHistMap;
 std::map<std::shared_ptr<angleHistData>, std::vector<TH1F*>, pawian::Collection::SharedPtrLess > _angleGJFitHistMap;

 virtual void initRootStuff()=0;

private:
  float _weightToWrite;
  std::vector<Particle*> _fsParticles;
  std::map<std::string, std::shared_ptr<TLorentzVector> > _fourVecMap;

};

