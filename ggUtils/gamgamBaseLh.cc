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

// gamgamBaseLh class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "ggUtils/gamgamBaseLh.hh"
#include "ggUtils/ggReaction.hh"
#include "ggUtils/GGChannelEnv.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "PwaUtils/FsParticleProjections.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"
#include "ConfigParser/ParserBase.hh"
#include "ConfigParser/ggParser.hh"
#include "Utils/FunctionUtils.hh"
#include "Utils/IdStringMapRegistry.hh"
#include "qft++Extension/PawianUtils.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

gamgamBaseLh::gamgamBaseLh(ChannelID channelID) :
  AbsLh(channelID)
  ,_ggChannelEnv(std::static_pointer_cast<GGChannelEnv> (GlobalEnv::instance()->GGChannel(channelID)))
  ,_useProdDynamics(false)
  ,_withDecMom(false)
  ,_sProdExp(0.0)
{
  initialize();
}


gamgamBaseLh::~gamgamBaseLh()
{;
}


complex<double> gamgamBaseLh::calcSpinDensity(Spin M1, Spin M2, std::string& nameDec, EvtData* theData){

  complex<double> result(0.,0.);

  return result;
}




complex<double> gamgamBaseLh::calcProdPartAmp(Spin lamX, Spin lamDec, std::string nameDec, EvtData* theData,
					     std::map <std::shared_ptr<const JPCLS>,
					     std::vector< std::shared_ptr<AbsXdecAmp> >,
					     pawian::Collection::SharedPtrLess > pbarpAmps){
   complex<double> resultAmp(0.,0.);

   return resultAmp;
}


double gamgamBaseLh::calcEvtIntensity( EvtData* theData, std::shared_ptr<AbsPawianParameters> fitPar){

  double result=0.;

  std::vector< std::shared_ptr<AbsXdecAmp> >::iterator itDecAll;
  for (itDecAll=_decAmps.begin(); itDecAll!=_decAmps.end(); ++itDecAll){
    (*itDecAll)->calcDynamics(theData);
  }

  const std::vector< std::vector<Spin> >& spinProjections=_fsParticleProjections->spinProjections();
  
  for (unsigned int projId=0; projId<spinProjections.size(); ++projId){
    for (itDecAll=_decAmps.begin(); itDecAll!=_decAmps.end(); ++itDecAll){
      //(*itDecAll)->setSpinProjections(projId);
      const std::vector<Spin>& currentSpinProjection=spinProjections.at(projId);
      (*itDecAll)->setSpinProjections(currentSpinProjection, projId);
    }
    
    complex<double> lamX0Amp(0.,0.);
    
    std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
    Spin lamX=0;
    for( itDec=_decAmps.begin(); itDec!=_decAmps.end(); ++itDec){
      complex<double> currentDecAmp=(*itDec)->XdecAmp(lamX, theData);
      lamX0Amp+=currentDecAmp;
    }


    complex<double> lamXp2Amp(0.,0.);
    lamX=2;
    for( itDec=_decAmps.begin(); itDec!=_decAmps.end(); ++itDec){
      complex<double> currentDecAmp=(*itDec)->XdecAmp(lamX, theData);
      lamXp2Amp+=currentDecAmp;
    }

    complex<double> lamXm2Amp(0.,0.);
    lamX=-2;
    for( itDec=_decAmps.begin(); itDec!=_decAmps.end(); ++itDec){
      complex<double> currentDecAmp=(*itDec)->XdecAmp(lamX, theData);
      lamXm2Amp+=currentDecAmp;
    }

   
    if(_useCohPhasespace){
      lamX0Amp += std::polar( fitPar->Value(_CohPhasespaceKey+"Mag"), 
			      fitPar->Value(_CohPhasespaceKey+"Phi") );
      lamXp2Amp += std::polar( fitPar->Value(_CohPhasespaceKey+"Mag"), 
			       fitPar->Value(_CohPhasespaceKey+"Phi") );
      lamXm2Amp += std::polar( fitPar->Value(_CohPhasespaceKey+"Mag"), 
			       fitPar->Value(_CohPhasespaceKey+"Phi") );
    }
 
    result += norm(lamX0Amp) + norm(lamXp2Amp) + norm(lamXm2Amp);
  }
  
  if(_usePhasespace) result+=fitPar->Value(_phasespaceKey);
  if(_useProdDynamics){
    result*=fabs(_dyn->eval(theData,0).real());
    double currentMass=_dyn->currentMass(theData);
    result*=pow(currentMass*currentMass,_sProdExp);    
    if(_withDecMom){
      double currentMass=_dyn->currentMass(theData);
      double currentBreakupMom=(PawianQFT::breakupMomQDefault(currentMass, _decMass1, _decMass2)).real();
      result*=currentBreakupMom;
    }
  }
  result *= fitPar->Value(_channelScaleParam);
  return result;
}

void gamgamBaseLh::print(std::ostream& os) const{

}

void  gamgamBaseLh::initialize(){
  std::vector< std::shared_ptr<AbsDecay> > theDecs = _ggChannelEnv->prodDecayList()->getList();
  std::vector< std::shared_ptr<AbsDecay> >::iterator it;
  for (it=theDecs.begin(); it!=theDecs.end(); ++it){
    InfoMsg << "theDecs->name: " << (*it)->name() << endmsg;
    std::shared_ptr<AbsXdecAmp> currentAmp=XdecAmpRegistry::instance()->getXdecAmp(_channelID, 
										   (*it)->absDecPtr());
    _decAmps.push_back(currentAmp);
  }

 std::string dynString=_ggChannelEnv->parser()->productionDynamics();
  if(dynString=="woFormPol" || dynString=="FormPol0" || dynString=="FormPol1" || dynString=="FormPol2"){
    _useProdDynamics=true;
    std::shared_ptr<FormationDecay> motherFormDec = _ggChannelEnv->reaction()->motherProdDec();
    std::vector<std::string> additionalStringVecDummy;
    motherFormDec->enableDynamics(dynString, additionalStringVecDummy);
   _dyn=motherFormDec->getDynamics();
  }
  
  ggParser* theggParser = dynamic_cast<ggParser*>(_ggChannelEnv->parser());
  std::string decMomStr=theggParser->decMomParamsStr();
  std::stringstream decMomParamsStringStr;
  decMomParamsStringStr << decMomStr;
  
  std::string withDecMomStr;
  decMomParamsStringStr >> withDecMomStr;
  int withDecMom=atof(withDecMomStr.c_str());
  if (withDecMom!=0) _withDecMom=true;

  std::string decMass1Str;
  decMomParamsStringStr >> decMass1Str;
  _decMass1=atof(decMass1Str.c_str());

  std::string decMass2Str;
  decMomParamsStringStr >> decMass2Str;
  //  _decMass1=atof(decMass2Str.c_str());
  _decMass2=atof(decMass2Str.c_str()); 

  _sProdExp = theggParser->sProdExponent();
}

void gamgamBaseLh::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  if(_useProdDynamics) _dyn->fillDefaultParams(fitPar);
 AbsLh::fillDefaultParams(fitPar);
}

void gamgamBaseLh::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  if(_useProdDynamics) _dyn->updateFitParams(fitPar);
  AbsLh::updateFitParams(fitPar);
}

