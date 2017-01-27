//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// pipiScatteringBaseLh class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "pipiScatteringUtils/pipiScatteringBaseLh.hh"
#include "pipiScatteringUtils/PiPiScatteringChannelEnv.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "ErrLogger/ErrLogger.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


pipiScatteringBaseLh::pipiScatteringBaseLh(ChannelID channelID) :
  AbsLh(channelID)
{
  initialize();
}


pipiScatteringBaseLh::~pipiScatteringBaseLh()
{;
}


complex<double> pipiScatteringBaseLh::calcSpinDensity(Spin M1, Spin M2, std::string& nameDec, EvtData* theData){

  complex<double> result(0.,0.);

  return result;
}




complex<double> pipiScatteringBaseLh::calcProdPartAmp(Spin lamX, Spin lamDec, std::string nameDec, EvtData* theData,
					     std::map <std::shared_ptr<const JPCLS>,
					     std::vector< std::shared_ptr<AbsXdecAmp> >,
					     pawian::Collection::SharedPtrLess > resAmps){
   complex<double> resultAmp(0.,0.);

   return resultAmp;
}


double pipiScatteringBaseLh::calcEvtIntensity(EvtData* theData, std::shared_ptr<AbsPawianParameters> fitPar){
  // chi2 fit must be included here
  // each point represtes one event???
  double result=0.;

  return result;

}

void pipiScatteringBaseLh::print(std::ostream& os) const{

}


void  pipiScatteringBaseLh::initialize(){
  const std::shared_ptr<PiPiScatteringChannelEnv> pipiScatteringEnv=std::static_pointer_cast<PiPiScatteringChannelEnv>(GlobalEnv::instance()->PiPiScatteringChannel(_channelID));

  std::vector<std::shared_ptr<AbsDecay> > absDecList=pipiScatteringEnv->absDecayList()->getList();
  _XdecAmp=XdecAmpRegistry::instance()->getXdecAmp(_channelID, absDecList.at(0)->absDecPtr());
  _decAmps.push_back(_XdecAmp);
}


void pipiScatteringBaseLh::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  _XdecAmp->fillDefaultParams(fitPar);
}

void pipiScatteringBaseLh::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  _XdecAmp->updateFitParams(fitPar);
}


