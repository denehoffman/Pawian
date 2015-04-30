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

#include <fstream>

#include "FitParams/AbsPawianParameters.hh"
#include "ErrLogger/ErrLogger.hh"




AbsPawianParameters::AbsPawianParameters()
{
}



AbsPawianParameters::~AbsPawianParameters(){;}

void AbsPawianParameters::print(std::ostream& os){
  std::vector<std::string> nameVec=ParamNames();

  std::vector<std::string>::const_iterator it;
  for(it=nameVec.begin(); it!=nameVec.end() ; ++it){
    //    std::cout << *it;
    // std::cout << "\tidx: " << Index(*it) << std::endl;
    int idx=Index(*it);
    os << GetName(idx) << "\t" << Value(idx) << "\t" << Error(idx);
    if(HasLimits(idx)) os << "\t" << LowerLimit(idx) << "\t" << UpperLimit(idx);
    os << "\n"; 
    }
}

void AbsPawianParameters::SetAllValues(const std::vector<double>& values){
  if(values.size() != Params().size()){
    Alert << "different number of parameter!!!"
	  << "\nvalues.size(): " << values.size() << "\tParams().size(): " << Params().size() <<endmsg;
    exit(0);
  }

  for(unsigned int i=0; i<values.size(); ++i){
    SetValue(i, values.at(i));
  }
}


