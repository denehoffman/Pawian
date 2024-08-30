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

#include <math.h>
#include <stdio.h>
#include <iomanip>
#include "Minuit2/FCNBase.h"
#include "Minuit2/FCNGradientBase.h"

#include "MinFunctions/PwaFcnServer.hh"
#include "FitParams/ParamDepHandler.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/NetworkServer.hh"
#include "Utils/PawianConstants.hh"
#include "ConfigParser/ParserBase.hh"
#include "ErrLogger/ErrLogger.hh"

using namespace ROOT::Minuit2;

template<typename T>
PwaFcnServer<T>::PwaFcnServer(std::shared_ptr<NetworkServer> netServer) :
  AbsFcn<T>()
  , _networkServerPtr(netServer)
  , _numStepSize(std::sqrt(std::numeric_limits<double>::epsilon()*150.))
{
  this->_defaultPawianParms = GlobalEnv::instance()->defaultPawianParams();
  this->_currentPawianParms = GlobalEnv::instance()->startPawianParams();
}

template<typename T>
PwaFcnServer<T>::~PwaFcnServer()
{
}

template<typename T>
double PwaFcnServer<T>::operator()(const std::vector<double>& par) const
{
  //double result=0.;
  this->_currentPawianParms->SetAllValues(par);
  ParamDepHandler::instance()->ApplyDependencies(this->_currentPawianParms);

  return collectLH(); 
}

template<typename T>
double PwaFcnServer<T>::collectLH() const{
  double result=0;
  std::map<ChannelID, LHData> theLHDataMap;
  std::ostringstream output;
  std::ostringstream outputLHDump;
  bool lhPrint=false;
  if(this->_fcnCounter%GlobalEnv::instance()->parser()->stepSizeLhPrint() == 0) lhPrint=true;
  bool lhDump=false;
  if(this->_fcnCounter%GlobalEnv::instance()->parser()->stepSizeLhDump() == 0) lhDump=true;
  
  _networkServerPtr->BroadcastParams(this->_currentPawianParms->Params());
  if(!_networkServerPtr->WaitForLH(theLHDataMap))
    result = 0;
  else{
    if(lhDump && theLHDataMap.size() > 1){
      outputLHDump << result << "\t";
    }
      // Add LLHs of different channels
    if(lhPrint) output << "current LH = ";
      for(auto it = theLHDataMap.begin(); it!=theLHDataMap.end();++it){
         (*it).second.weightSum = _networkServerPtr->weightSum((*it).first);
	 (*it).second.squaredWeightSum = _networkServerPtr->squaredWeightSum((*it).first);
         (*it).second.num_mc = _networkServerPtr->numMCs((*it).first);
         double channelLH = AbsLh::mergeLogLhData((*it).second, (*it).first);
         result += channelLH;
	 if(lhPrint) output << std::setprecision(16) << channelLH << "\t";
	 if(lhDump) outputLHDump << std::setprecision(16) << channelLH << "\t";
      }
      if(lhPrint && theLHDataMap.size() > 1){
         output << "sum = " << result;
      }
  }

  if(lhPrint){
    InfoMsg << output.str() << endmsg;
  }
  
  if(this->_fcnCounter%GlobalEnv::instance()->parser()->stepSizeTimer() == 0) this->printTimer();
  if(this->_fcnCounter%GlobalEnv::instance()->parser()->stepSizeParamsPrint() == 0) this->printFitParams(this->_currentPawianParms);
  if(this->_fcnCounter%GlobalEnv::instance()->parser()->stepSizeParamsDump() == 0) this->dumpFitParams(this->_currentPawianParms);
  if(lhDump) this->dumpLhVals(outputLHDump.str());
  this->_fcnCounter++;
  return result;
}

template<typename T>
std::vector<double> PwaFcnServer<T>::Gradient(const std::vector<double>& par) const{
  std::vector<double> resultVec;
  resultVec.resize(par.size());
  this->_currentPawianParms->SetAllValues(par);
  ParamDepHandler::instance()->ApplyDependencies(this->_currentPawianParms);
  double LHBase=collectLH();
  double epsilon;

  for(unsigned int i=0; i<par.size(); ++i){
    if(this->_currentPawianParms->IsFixed(i)) resultVec.at(i)=0.;
    else{
      double currentVal=this->_currentPawianParms->Value(i);
      epsilon=_numStepSize*std::abs(currentVal);
      if ((this->_currentPawianParms->GetName(i)).substr( (this->_currentPawianParms->GetName(i)).length() - 3 ) == "Phi"){
	epsilon=_numStepSize*PawianConstants::pi;
      }
      else if (std::abs(currentVal)<1.e-10)  epsilon=_numStepSize*1.e-10;
      double dx=(currentVal+epsilon)-currentVal;
      if(this->_currentPawianParms->HasLimits(i) && std::abs(currentVal-this->_currentPawianParms->UpperLimit(i))<1.e-6){
	this->_currentPawianParms->SetValue(i, currentVal-epsilon);
        ParamDepHandler::instance()->ApplyDependencies(this->_currentPawianParms);
        double currentLH=collectLH();
	resultVec.at(i)=(LHBase-currentLH)/dx;
	//           InfoMsg << "resultVecLow.at(" << i << ")= " << resultVec.at(i) << endmsg;
      }
      else{
	this->_currentPawianParms->SetValue(i, currentVal+epsilon);
      	ParamDepHandler::instance()->ApplyDependencies(this->_currentPawianParms);
       	double currentLH=collectLH();
	resultVec.at(i)=(currentLH-LHBase)/dx;
	//           InfoMsg << "resultVecHigh.at(" << i << ")= " << resultVec.at(i) << endmsg;
      }
      this->_currentPawianParms->SetValue(i, currentVal);
    }
  }
  return resultVec;  
}

template PwaFcnServer<FCNBase>::PwaFcnServer(std::shared_ptr<NetworkServer>);
template PwaFcnServer<FCNBase>::~PwaFcnServer();
template PwaFcnServer<FCNGradientBase>::PwaFcnServer(std::shared_ptr<NetworkServer>);
template PwaFcnServer<FCNGradientBase>::~PwaFcnServer();




