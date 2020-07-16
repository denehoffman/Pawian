//************************************************************************//
//									  //
//  Copyright 2017 Markus Kuhlmann (mkuhlmann@ep1.rub.de)	          //
//  	      	   - Ruhr-Universität Bochum 				  //
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
#include <memory>

#include "Minuit2/FCNBase.h"

class BBUnstableParFit;

namespace ROOT {
namespace Minuit2 {

class BBUnstableParFitFcn : public FCNBase {

public:

  BBUnstableParFitFcn(std::shared_ptr<BBUnstableParFit>);
  virtual ~BBUnstableParFitFcn();

  double operator()(const std::vector<double>& par) const;

  double Up() const;

private:
  std::shared_ptr<BBUnstableParFit> _minFitPtr;
};

}  // namespace Minuit2
}  // namespace ROOT
