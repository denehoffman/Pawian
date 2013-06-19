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
#include <memory>

#include "Minuit2/FCNBase.h"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "Minuit2/MnUserParameters.h"

namespace ROOT {
  namespace Minuit2 {
    class AbsFcn : public FCNBase {

    public:
      AbsFcn(std::shared_ptr<FitParamsBase> fitParamsBase, std::string suffix="");
      virtual ~AbsFcn();

      virtual double operator()(const std::vector<double>& par) const=0;
      virtual double Up() const;


    protected:
      std::shared_ptr<FitParamsBase> _fitParamsBasePtr;
      mutable unsigned int _fcnCounter;
       fitParams _defaultFitValParms;
      fitParams _defaultFitErrParms;
      std::string _currentResFileName;

      virtual void printTimer() const;
      virtual void printFitParams(const std::vector<double>& par) const;
      virtual void dumpFitParams(const std::vector<double>& par) const;

    private:
    };
  }  // namespace Minuit2
}  // namespace ROOT
