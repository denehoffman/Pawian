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
#include "FitParams/ParamFactory.hh"
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

std::shared_ptr<AbsPawianParameters> AbsPawianParameters::paramsWithSameOrder(std::shared_ptr<AbsPawianParameters> toSort){
  if(Params().size() != toSort->Params().size() ){
    Alert << "parameter lists exhibit differnet sizes: Params().size(): " << Params().size() 
	  << " toSort->Params().size(): " << toSort->Params().size() << endmsg;
    exit(1); 
  }
  
  std::shared_ptr<AbsPawianParameters> sortedList=ParamFactory::instance()->getParametersPointer(toSort->type());

  for(unsigned int id=0; id<Params().size(); ++id){
    const std::string currentName=GetName(id);
    unsigned int toSortId=toSort->Index(currentName);
    if(toSortId >= Params().size()){
      Alert << "toSort list does not contain parameter with the name " << currentName << endmsg;
      exit(1); 
    }
    sortedList->Add(currentName, toSort->Value(toSortId), toSort->Error(toSortId));
    if(toSort->HasLimits(toSortId)){
      sortedList->SetLimits(id, toSort->LowerLimit(toSortId), toSort->UpperLimit(toSortId));
    }
    if(toSort->IsFixed(toSortId)) sortedList->Fix(id);    
  }

  return sortedList;
}


