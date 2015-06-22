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
#include <iomanip>

#include "FitParams/AbsPawianParameters.hh"
#include "FitParams/ParamFactory.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Minuit2/MnUserParameters.h"



AbsPawianParameters::AbsPawianParameters()
{
}



AbsPawianParameters::~AbsPawianParameters(){;}

std::shared_ptr<MnUserParameters> AbsPawianParameters::mnUserParametersPtr(){
  std::shared_ptr<MnUserParameters> result(new MnUserParameters());
  for(unsigned int i=0; i<Params().size(); ++i){
    const std::string currentName=GetName(i);
    const double currentVal=Value(i);
    const double currentErr=Error(i);
    if(HasLimits(i)){
      const double currentLoLimit=LowerLimit(i);
      const double currentUpLimit=UpperLimit(i);
      result->Add(currentName, currentVal, currentErr, currentLoLimit, currentUpLimit);  
    }
    else{
      result->Add(currentName, currentVal, currentErr);  
    }

    if(IsFixed(i)) result->Fix(currentName); 
  }
  return result;
}

void AbsPawianParameters::print(std::ostream& os, bool extended){
  std::vector<std::string> nameVec=ParamNames();

  std::vector<std::string>::const_iterator it;
  if(!extended){
    for(it=nameVec.begin(); it!=nameVec.end() ; ++it){
      //    std::cout << *it;
      // std::cout << "\tidx: " << Index(*it) << std::endl;
      int idx=Index(*it);
      os << GetName(idx) << "\t" << Value(idx) << "\t" << Error(idx);
      if(HasLimits(idx)) os << "\t" << LowerLimit(idx) << "\t" << UpperLimit(idx);
      os << "\n"; 
    }
  }
  else{
    unsigned int maxNameLength=10;
    for(it=nameVec.begin(); it!=nameVec.end() ; ++it){
      if ( it->size()>maxNameLength ) maxNameLength=it->size();
    }
    
    int smallwidth = 5;
    int mediumwidth = 15;
    //  int largewidth = 25;
    os << std::setw(smallwidth) << "no " << std::setw(maxNameLength+1) << "Name" << std::setw(mediumwidth) << "Value" << std::setw(smallwidth+1) << "type" << std::setw(mediumwidth) << "Error +-" << std::setw(mediumwidth) << "Low" << std::setw(mediumwidth) << "High\n";
    for(it=nameVec.begin(); it!=nameVec.end() ; ++it){
      int idx=Index(*it);
      os << std::setw(smallwidth-1) << idx << " " << std::setw(maxNameLength+1) << GetName(idx) << std::setw(mediumwidth) << Value(idx);
      if (IsFixed(idx)) os << std::setw(smallwidth+1) << "fixed";
      else{
	os << std::setw(smallwidth+1) << "free" << std::setw(mediumwidth) << Error(idx);
	if(HasLimits(idx)) os << std::setw(mediumwidth) << LowerLimit(idx) << std::setw(mediumwidth) << UpperLimit(idx);
      }
      os << "\n";
    }
  }
}

void AbsPawianParameters::SetAllValues(const std::vector<double>& values){
  if(values.size() != Params().size()){
    Alert << "different number of parameter!!!"
	  << "\nvalues.size(): " << values.size() << "\tParams().size(): " << Params().size() <<endmsg;
    exit(0);
  }

  for(unsigned int i=0; i<values.size(); ++i){
    if(values.at(i) != values.at(i)) {
       Alert << "Parameter vector returned from Minuit conains NAN! Exiting..." << endmsg;
       exit(0);
    }
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


