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

// gammapStates class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#include "gammapUtils/gammapStates.hh"
#include "Utils/MathUtils.hh"
#include "ErrLogger/ErrLogger.hh"

gammapStates::gammapStates():
  AbsStates(),
  _lmax(10),
  _pJPC(0.5, 1),
  _gammaJPC(1,-1,-1)
{
   calcStates();
}


gammapStates::gammapStates(int lmax):
  AbsStates(),
  _lmax(lmax),
  _pJPC(0.5, 1),
  _gammaJPC(1,-1,-1)
{
    calcStates();
}

gammapStates::~gammapStates(){
}

bool gammapStates::calcStates(){
  Spin jmax=_lmax+Spin(3./2);
  std::cout << "jmax: " << jmax << std::endl;
  std::cout << "lmax: " << _lmax << std::endl; 
 
  for (Spin j=Spin(1./2); j<=jmax; ++j){
    for (int p=-1; p<=1; p+=2){

      std::shared_ptr<const jpcRes> jpcPtr(new jpcRes(j,p,0));
      std::shared_ptr<const IGJPC> i12gjpcPtr(new IGJPC(jpcPtr, Spin(1./2), 0));
      std::shared_ptr<const IGJPC> i32gjpcPtr(new IGJPC(jpcPtr, Spin(3./2), 0));
      
      //now fill multupoles
      Spin j;
      Spin l;
      if(jpcPtr->P == pow(-1.,(int)(jpcPtr->J + 1./2))){ // 1/2-,3/2+,5/2-,...
	l=jpcPtr->J-Spin(1./2);
	if(l<=_lmax){
	  fillVec(jpcPtr, _alljpcRes);
	  fillVec(i12gjpcPtr, _allIGjpcRes);
	  fillVec(i32gjpcPtr, _allIGjpcRes);
	  if(jpcPtr->J < 3/2.){
	    j=1;
	    std::shared_ptr<const JPCLS> tmpJPClj( new JPCLS(jpcPtr,l, j) );
	    fillVec(tmpJPClj, _JPClj_States); //E1 multipole
	    fillVec(tmpJPClj,_elMpolJPClj_States);
	  }
	  else{
	    //first fill electric multipoles
	    j=jpcPtr->J+Spin(1./2);
	    std::shared_ptr<const JPCLS> tmpJPClj_el( new JPCLS(jpcPtr,l, j) );
	    fillVec(tmpJPClj_el, _JPClj_States); //Ej multipole 
	    fillVec(tmpJPClj_el,_elMpolJPClj_States);
	    
	    //now fill magnetic multipoles
	    j=jpcPtr->J-Spin(1./2);
	    std::shared_ptr<const JPCLS> tmpJPClj_mag( new JPCLS(jpcPtr,l, j) );
	    fillVec(tmpJPClj_mag, _JPClj_States); //Mj multipole
	    fillVec(tmpJPClj_mag, _magMpolJPClj_States);                   
	  }
	}
      }
      else{ // 1/2+,3/2-,5/2+,...
	l=jpcPtr->J+Spin(1./2);
	if(l<=_lmax){
	  fillVec(jpcPtr, _alljpcRes);
	  fillVec(i12gjpcPtr, _allIGjpcRes);
	  fillVec(i32gjpcPtr, _allIGjpcRes);
	  if(jpcPtr->J < 3/2.){
	    j=1;
	    std::shared_ptr<const JPCLS> tmpJPClj( new JPCLS(jpcPtr,l, j) );
	    fillVec(tmpJPClj, _JPClj_States); //M1 multipole
	    fillVec(tmpJPClj, _magMpolJPClj_States);  
	  }
	  else{
	    //first fill electric multipoles
	    j=jpcPtr->J-Spin(1./2);
	    std::shared_ptr<const JPCLS> tmpJPClj_el( new JPCLS(jpcPtr,l, j) );
	    fillVec(tmpJPClj_el, _JPClj_States); //Ej multipole
	    fillVec(tmpJPClj_el, _elMpolJPClj_States); 
	    
	    //now fill magnetic multipoles
	    j=jpcPtr->J+Spin(1./2);
	    std::shared_ptr<const JPCLS> tmpJPClj_mag( new JPCLS(jpcPtr,l, j) );
	    fillVec(tmpJPClj_mag, _JPClj_States); //Mj multipole
	    fillVec(tmpJPClj_mag, _magMpolJPClj_States);                     
	  }
	}
      }
    }
  }
  
  return true;
}


void gammapStates::print(std::ostream& os) const{
  os << "initital states of the p gamma reaction for Lmax = " << _lmax << " are: " << std::endl; 
  //  AbsStates::print(os);
   std::vector< std::shared_ptr<const JPCLS > >::const_iterator itJPCLS;
   std::vector< std::shared_ptr<const jpcRes > >::const_iterator itjpcRes;
   std::vector< std::shared_ptr<const IGJPC > >::const_iterator itIGJPC;
  os << "\n******** all JPC states ************ " << std::endl;
  for ( itjpcRes=_alljpcRes.begin(); itjpcRes!=_alljpcRes.end(); ++itjpcRes){
    (*itjpcRes)->print(os);
    os << std::endl;
  }

  os << "\n******** all IGJPC states ************ " << std::endl;
  for ( itIGJPC=_allIGjpcRes.begin(); itIGJPC!=_allIGjpcRes.end(); ++itIGJPC){ 
    (*itIGJPC)->print(os);
    os << std::endl;
  }

 os << "\n******** multipoles JP states  with l j(gamma)************ " << std::endl;
 for ( itJPCLS=_JPClj_States.begin(); itJPCLS!=_JPClj_States.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << std::endl;
  }

 os << "\n******** electric multipoles JP states  with l j(gamma)************ " << std::endl;
 for ( itJPCLS=_elMpolJPClj_States.begin(); itJPCLS!=_elMpolJPClj_States.end(); ++itJPCLS){
   os << "\nE" << (*itJPCLS)->S << "\t";
    (*itJPCLS)->print(os);
    os << std::endl;
  }

 os << "\n******** magnetic multipoles JP states  with l j(gamma)************ " << std::endl;
 for ( itJPCLS=_magMpolJPClj_States.begin(); itJPCLS!=_magMpolJPClj_States.end(); ++itJPCLS){
   os << "\nM" << (*itJPCLS)->S << "\t";
    (*itJPCLS)->print(os);
    os << std::endl;
  }

}

