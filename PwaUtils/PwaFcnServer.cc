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
#include <boost/timer/timer.hpp>

#include "Minuit2/MnUserParameters.h"

#include "PwaUtils/PwaFcnServer.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/NetworkServer.hh"
#include "ErrLogger/ErrLogger.hh"

using namespace ROOT::Minuit2;

boost::timer::cpu_timer theTimer;

PwaFcnServer::PwaFcnServer(boost::shared_ptr<AbsLh> absLh, boost::shared_ptr<FitParamsBase> fitParamsBase, boost::shared_ptr<NetworkServer> netServer, std::string suffix) :
  _absLhPtr(absLh)
  , _fitParamsBasePtr(fitParamsBase)
  , _networkServerPtr(netServer)
  , _fcnCounter(0)
  , _currentResFileName("currentResult"+suffix+".dat")
{
   if (0==_absLhPtr) { Alert << "AbsLh* _absLhPtr pointer is 0 !!!!" << endmsg; exit(1); }
   _absLhPtr->getDefaultParams(_defaultFitValParms, _defaultFitErrParms);
  
}

PwaFcnServer::~PwaFcnServer()
{
}

double PwaFcnServer::operator()(const std::vector<double>& par) const
{
  double result=0;
  
  LHData theLHData;
  _networkServerPtr->BroadcastParams(par);
  if(!_networkServerPtr->WaitForLH(theLHData.logLH_data, theLHData.weightSum, theLHData.LH_mc))
    result = 0;
  else
    result = _absLhPtr->mergeLogLhData(theLHData, _networkServerPtr->numMCs());
  
  
  _fcnCounter++;
  
  if(_fcnCounter%20 == 0){
    theTimer.stop();
    boost::timer::cpu_times elapsed(theTimer.elapsed());
    if(elapsed.wall > 0){
      Info << "Wall time: " << elapsed.wall / 1E9 << "s User: "
	   << elapsed.user/1E9 << "s System: " << elapsed.system/1E9 << "s\n" << endmsg;
    }
    theTimer.start();
  }
  
  if (  _fcnCounter%100 == 0) {
    fitParams theFitParmValTmp=_defaultFitValParms;
    _fitParamsBasePtr->getFitParamVal(par, theFitParmValTmp);
    _fitParamsBasePtr->printParams(theFitParmValTmp);
    if (  _fcnCounter%200 == 0) {
      std::ofstream theStream (_currentResFileName.c_str());
      _fitParamsBasePtr->dumpParams(theStream, theFitParmValTmp, (fitParams&)_defaultFitErrParms);
    }
  }
  
  return result;
}

double PwaFcnServer::Up() const 
{
return .5;
}



