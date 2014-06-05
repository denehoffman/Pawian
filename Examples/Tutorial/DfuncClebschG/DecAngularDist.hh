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

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cassert>

#include <memory>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

class jpcRes;
class TFile;
class TH1F;

class DecAngularDist {

public:

  // create/copy/destroy:

  ///Constructor 
  DecAngularDist(std::shared_ptr<const jpcRes> motherJPC, std::shared_ptr<const jpcRes> daughter1JPC, std::shared_ptr<const jpcRes> daughter2JPC, bool sameDaughters=false);


  /** Destructor */
  virtual ~DecAngularDist();

  // Getters:

protected:


private:
  std::shared_ptr<const jpcRes> _motherJPC;
  std::shared_ptr<const jpcRes> _daughter1JPC;
  std::shared_ptr<const jpcRes> _daughter2JPC;
  std::vector<LS> _allLSs;
  std::vector<LS> _allowedCParityLSs;
  std::vector<LS> _forbiddenCParityLSs;
  bool _sameDaughters;
  LeviCivitaTensor _lctTensor;
  MetricTensor _metricTensor;
  PolVector _polMother;
  PolVector _polDaughter1;
  PolVector _polDaughter2;
  double _qCMSReal;
  Vector4<double> _mother4Vec;

  TFile* _theTFile;

  std::map<Spin, std::map<Spin, std::map<Spin, TH1F*> > > _tensorJMM1M2AmpThetaHistMap;
  std::map<Spin, std::map<Spin, std::map<Spin, TH1F*> > > _canoJMM1M2AmpThetaHistMap; 

  std::map<Spin, TH1F*> _canoJMAmpThetaHistMap;
  TH1F* _canoJMAmpAllThetaHist;
  TH1F* _tensorJMAmpAllThetaHist;
  TH1F* _heliJMAmpAllThetaHist; 

  std::map<Spin, std::map<Spin, TH1F*> > _canoLSincohMAmpAllThetaHist;
  std::map<Spin, std::map<Spin, TH1F*> > _tensorLSincohMAmpAllThetaHist;

  std::map<Spin, std::map<Spin, std::map<Spin, std::map<Spin, std::map<Spin, TH1F*> > > > > _tensorJLSMM1M2AmpThetaHistMap;
  std::map<Spin, std::map<Spin, std::map<Spin, std::map<Spin, std::map<Spin, TH1F*> > > > > _canoJLSMM1M2AmpThetaHistMap;

  TH1F* _heliJAmpAllThetaHist;
  TH1F* _heliLSJAmpAllThetaHist;
  TH1F* _canoJAmpAllThetaHist;
  TH1F* _tensorJAmpAllThetaHist;
  TH1F* _tensorJMM1M2AmpAllThetaHist;
  TH1F* _canoJMM1M2AmpAllThetaHist;
  
  void calcChi12(Tensor<complex<double> >& s12SpinProjector, Tensor<complex<double> >& epsilonDaughter1Project, Tensor<complex<double> >& epsilonDaughter2Project, Vector4<double>& mother_4Vec, Tensor<complex<double> >& result);

  void calcLSpart(OrbitalTensor& orbTensor, Tensor<complex<double> >& chi12, Spin spin12, Vector4<double>& mother_4Vec, Tensor<complex<double> >& result);

  void fillTensorAmps(double cosT);
  void fillCanoAmps(double cosT);
  void fillCanoAmps1(double cosT);
  void fillheliAmps(double cosT);
};

