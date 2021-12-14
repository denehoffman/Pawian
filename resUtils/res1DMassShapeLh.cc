//************************************************************************//
//									  //
//  Copyright 2021 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// res1DMassShapeLh class definition file. -*- C++ -*-
// Copyright 2021 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "resUtils/res1DMassShapeLh.hh"
#include "resUtils/resReaction.hh"
#include "resUtils/ResChannelEnv.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


res1DMassShapeLh::res1DMassShapeLh(ChannelID channelID) :
  resBaseLh(channelID)
{

  std::vector< std::shared_ptr<AbsXdecAmp> >::iterator itDecAll;
  for (itDecAll=_decAmps.begin(); itDecAll!=_decAmps.end(); ++itDecAll){
    //InfoMsg << "amp name: " << (*itDecAll)->name() << endmsg;
    std::shared_ptr<AbsDynamics> grandMaDyn=(*itDecAll)->getDyn();
    std::shared_ptr<AbsXdecAmp> xDecDaughter1=(*itDecAll)->absXDecAmpDaughter1();
    if(xDecDaughter1) {
      std::shared_ptr<AbsDynamics> currentDyn=xDecDaughter1->getDyn();
      currentDyn->setGrandMaDyn(grandMaDyn);
      _dynVecs.push_back(currentDyn);
    }
    std::shared_ptr<AbsXdecAmp> xDecDaughter2=(*itDecAll)->absXDecAmpDaughter2();
    if(xDecDaughter2) {
      std::shared_ptr<AbsDynamics> currentDyn=xDecDaughter2->getDyn();
      currentDyn->setGrandMaDyn(grandMaDyn);
      _dynVecs.push_back(currentDyn);
    }
  }
  
  if (_dynVecs.size()>1){
    Alert << "mode with prodFormalism=1DMassFit supports only one decay dynamics" << endmsg;
    exit(1);}
  
}


res1DMassShapeLh::~res1DMassShapeLh()
{;
}


double res1DMassShapeLh::calcEvtIntensity(EvtData* theData, std::shared_ptr<AbsPawianParameters> fitPar){
  AbsXdecAmp* dummyAmp=0;
  complex<double> result(0.,0.);

  std::vector< std::shared_ptr<AbsDynamics> >::iterator itDyn;
  for (itDyn=_dynVecs.begin(); itDyn!=_dynVecs.end(); ++itDyn){
    result += (*itDyn)->eval(theData, dummyAmp);
  }

  result *= fitPar->Value(_channelScaleParam);
  return norm(result);

}

void res1DMassShapeLh::print(std::ostream& os) const{

}


void res1DMassShapeLh::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  
  fitPar->Add(_channelScaleParam, 1., 0.01);
  fitPar->SetLimits(_channelScaleParam, 0., 20.);
  
  std::vector< std::shared_ptr<AbsDynamics> >::iterator itDyn;
  for (itDyn=_dynVecs.begin(); itDyn!=_dynVecs.end(); ++itDyn){
    (*itDyn)->fillDefaultParams(fitPar);
    InfoMsg << "dyn type: " << (*itDyn)->type() << endmsg;
    InfoMsg << "dyn name: " << (*itDyn)->name() << endmsg;
    std::vector<std::string> theParNames=(*itDyn)->paramNames();
    std::vector<std::string>::iterator itParName;
    for (itParName=theParNames.begin(); itParName!=theParNames.end(); ++itParName){
      InfoMsg << "param name: " << (*itParName) << endmsg;
    }  
  }

}

void res1DMassShapeLh::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){

 std::vector< std::shared_ptr<AbsDynamics> >::iterator itDyn;
  for (itDyn=_dynVecs.begin(); itDyn!=_dynVecs.end(); ++itDyn){
    (*itDyn)->updateFitParams(fitPar);
  }
}


