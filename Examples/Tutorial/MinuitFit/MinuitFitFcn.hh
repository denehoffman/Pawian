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

#ifndef _MinuitFitFcn_H
#define _MinuitFitFcn_H

#include <iostream>
#include <fstream>
//#include <string>
#include <vector>
#include <memory>
//#include <cassert>
#include "Minuit2/FCNBase.h"

class MinuitFit;

namespace ROOT {

   namespace Minuit2 {
class MinuitFitFcn : public FCNBase {

public:

  MinuitFitFcn(std::shared_ptr<MinuitFit>);
  virtual ~MinuitFitFcn();

  double operator()(const std::vector<double>& par) const;

  double Up() const;

private:
  std::shared_ptr<MinuitFit> _minFitPtr;
};
  }  // namespace Minuit2

}  // namespace ROOT


#endif 
