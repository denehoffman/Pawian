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
#include <vector>
#include <map>
#include <string>
#include <memory>

#include "FitParams/AbsFitParColStreamer.hh"
#include "FitParams/FitParColBase.hh"


class StreamFitParmsBase : public AbsFitParColStreamer {

public:
  StreamFitParmsBase(std::string&);
  virtual ~StreamFitParmsBase();

  fitParCol getFitParamVal() { return _paramVal;}
  fitParCol getFitParamErr() { return _paramErr;}

  virtual void fillParamMap() {return;}

protected:
  virtual void fillParams();

  virtual void fillJPCAmps(mapStrJPC& valMap, mapStrJPC& errMap,
			  const std::string& suffix);
 virtual void fillLSAmps(mapStrLS& valMap, mapStrLS& errMap,
			  const std::string& suffix);
  virtual void fillJPCLamLamAmps(mapStrJPCLamLam& valMap, mapStrJPCLamLam& errMap,
			      const std::string& suffix);
  virtual void fillJPCLSAmps(mapStrJPCLS& valMap, mapStrJPCLS& errMap,
			  const std::string& suffix);
  virtual void fillDoubles(mapStrDouble& valMap, mapStrDouble& errMap,
			   const std::string& suffix);
private:
  fitParCol _paramVal;
  fitParCol _paramErr;

  void fillParameter(std::map<int, double>& theValMap, std::map<int, double>& theErrMap,
		     std::string& suffix, int index);
};
