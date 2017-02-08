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
#include "MinFunctions/PwaFcnServer.hh"
#include "FitParams/ParamDepHandler.hh"

#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/NetworkServer.hh"
#include "ErrLogger/ErrLogger.hh"

using namespace ROOT::Minuit2;

PwaFcnServer::PwaFcnServer(std::shared_ptr<NetworkServer> netServer) :
  AbsFcn()
  , _networkServerPtr(netServer)
{
  _defaultPawianParms = GlobalEnv::instance()->defaultPawianParams();
   _currentPawianParms = GlobalEnv::instance()->startPawianParams();
}

PwaFcnServer::~PwaFcnServer()
{
}

double PwaFcnServer::operator()(const std::vector<double>& par) const
{
  double result=0;
  _currentPawianParms->SetAllValues(par);
  ParamDepHandler::instance()->ApplyDependencies(_currentPawianParms);

  std::map<ChannelID, LHData> theLHDataMap;
  _networkServerPtr->BroadcastParams(_currentPawianParms->Params());
  if(!_networkServerPtr->WaitForLH(theLHDataMap))
    result = 0;
  else{
      // Add LLHs of different channels
      std::ostringstream output;
      output << "current LH = ";
      for(auto it = theLHDataMap.begin(); it!=theLHDataMap.end();++it){
         (*it).second.weightSum = _networkServerPtr->weightSum((*it).first);
         (*it).second.num_mc = _networkServerPtr->numMCs((*it).first);
         double channelLH = AbsLh::mergeLogLhData((*it).second, (*it).first);
         result += channelLH;
	 output << std::setprecision(16) << channelLH << "\t";
      }
      if(theLHDataMap.size() > 1){
         output << "sum = " << result;
      }

      InfoMsg << output.str() << endmsg;
  }

  _fcnCounter++;
  if(_fcnCounter%20 == 0) printTimer();
  printFitParams(_currentPawianParms);

  dumpFitParams(_currentPawianParms);

  return result;
}




