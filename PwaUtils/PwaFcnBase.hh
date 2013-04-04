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
#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include "Minuit2/FCNBase.h"
//#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "Minuit2/MnUserParameters.h"

class AbsLh;
class NetworkServer;

namespace ROOT {
  namespace Minuit2 {
    class PwaFcnBase : public FCNBase {

    public:
      PwaFcnBase(boost::shared_ptr<AbsLh> absLh, 
		    boost::shared_ptr<FitParamsBase> fitParamsBase, std::string suffix="");
      virtual ~PwaFcnBase();

      double operator()(const std::vector<double>& par) const;
      double Up() const;
      void SetServerMode(boost::shared_ptr<NetworkServer> networkServerPtr);

    private:
      boost::shared_ptr<AbsLh> _absLhPtr;
      boost::shared_ptr<FitParamsBase> _fitParamsBasePtr;
      boost::shared_ptr<NetworkServer> _networkServerPtr;
      mutable unsigned int _fcnCounter;
      fitParams _defaultFitValParms;
      fitParams _defaultFitErrParms;
      std::string _currentResFileName;
      bool _serverMode;
    };
  }  // namespace Minuit2
}  // namespace ROOT
