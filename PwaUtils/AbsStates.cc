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

// AbsStates class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include "PwaUtils/AbsStates.hh"
#include "Utils/MathUtils.hh"

AbsStates::AbsStates(){
}


AbsStates::~AbsStates(){
}

// std::vector< boost::shared_ptr<const JPCLSM> > AbsStates::extractJPCLSMStates(std::vector< boost::shared_ptr<const jpcRes> >& theJPCStates) const{

//   std::vector< boost::shared_ptr<const JPCLSM> > result;

//   std::vector< boost::shared_ptr<const jpcRes> >::const_iterator itJPC;
//   std::vector< boost::shared_ptr<const JPCLSM> >::const_iterator itAllStates;

//   for ( itJPC=theJPCStates.begin(); itJPC!=theJPCStates.end(); ++itJPC){
//     const jpcRes* jpcRequest=(*itJPC).get();

//     for ( itAllStates=_allStates.begin(); itAllStates!=_allStates.end(); ++itAllStates){
//       const JPCLSM* jpcsmCurrent=(*itAllStates).get(); 
//       if (( *jpcsmCurrent) ==  (*jpcRequest)){
// 	result.push_back(*itAllStates);
// 	continue;
//       }
//     }
//   }
//   return result;
// }

// std::vector< boost::shared_ptr<const JPCSM> > AbsStates::extractJPCSMStates(std::vector< boost::shared_ptr<const jpcRes> >& theJPCStates) const{

//   std::vector< boost::shared_ptr<const JPCSM> > result;

//   std::vector< boost::shared_ptr<const jpcRes> >::const_iterator itJPC;
//   std::vector< boost::shared_ptr<const JPCSM> >::const_iterator itAllStates;

//   for ( itJPC=theJPCStates.begin(); itJPC!=theJPCStates.end(); ++itJPC){
//     const jpcRes* jpcRequest=(*itJPC).get();

//     for ( itAllStates=_allJPCSM.begin(); itAllStates!=_allJPCSM.end(); ++itAllStates){
//       const JPCSM* jpcsmCurrent=(*itAllStates).get(); 
//       if (( *jpcsmCurrent) ==  (*jpcRequest)){
// 	result.push_back(*itAllStates);
// 	continue;
//       }
//     }
//   }
//   return result;

// }

// std::vector< boost::shared_ptr<const JPCLS> > AbsStates::extractJPCLSStates(boost::shared_ptr<const jpcRes> theJPCState) const{
//   std::vector< boost::shared_ptr<const JPCLS> > result;

//   std::vector< boost::shared_ptr<const JPCLS> >::const_iterator it;
//   for ( it=_allJPCLS.begin(); it!=_allJPCLS.end(); ++it){
// //     const jpcRes* currentJPC= (jpcRes*) (*it).get(); 
// //     const jpcRes* currentJPC=dynamic_cast<const jpcRes*> ((*it).get());

// //     if ( *(currentJPC) == *(theJPCState.get())) result.push_back( (*it) );
//     if (*((*it).get())==*(theJPCState.get())) result.push_back( (*it) );
//   }
 
//   return result;
// }

void AbsStates::print(std::ostream& os) const{
  os << "AbsStates" << std::endl;
//   std::vector< boost::shared_ptr<const JPCLSM> >::const_iterator it;
//    for ( it=_allStates.begin(); it!=_allStates.end(); ++it){
//     if (0!= (*it)){
//       (*it)->print(os); 
// // 	os <<"\tL=" << (*it)->lsm.L <<"\tS=" << (*it)->lsm.S <<"\tlambda=" << (*it)->lsm.M
// // 	   <<"\tClebschGordan=" << (*it)->ClebschG 
// // 	   << std::endl;
      
//     }
//    }
}

