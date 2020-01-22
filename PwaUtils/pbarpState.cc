//************************************************************************//
//									  //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// pbarpState class definition file. -*- C++ -*-
// Copyright 2020 Bertram Kopf

#include "PwaUtils/pbarpState.hh"
#include "Utils/MathUtils.hh"
#include "ErrLogger/ErrLogger.hh"

pbarpState::pbarpState(std::shared_ptr<const jpcRes> theJPC):
  _jpc(theJPC),
  _pbarJPC(0.5, -1),
  _pJPC(0.5, 1),
  _ispbarpState(false),
  _isSinglet(false),
  _isTriplet(false),
  _isTriplet0(false), 
  _isTriplet1(false) 
{
  extractLSs();
}

pbarpState::~pbarpState(){
}

void pbarpState::extractLSs(){
  vector<LS> LScandidates=GetValidLS(_jpc->J, _jpc->P, _pbarJPC.J,_pbarJPC.P, _pJPC.J, _pJPC.P);
  int num_ls = (int) LScandidates.size();
  for(int ls = 0; ls < num_ls; ls++){
    bool accept=false;
    bool spinProjection0=false;
    bool spinProjection1=false; 
    Spin S= LScandidates.at(ls).S;
    Spin L= LScandidates.at(ls).L;
    int cparity(pow(-1,int(L+S)));
    if(cparity == _jpc->C){
      //check LS coupling
      for(Spin lampbar = -0.5; lampbar <= 0.5; ++lampbar){
	for(Spin lamp = -0.5; lamp <= 0.5; ++lamp){
	  double currentClebsch=Clebsch(L,0,S, lampbar-lamp, _jpc->J, lampbar-lamp);
          if( fabs(currentClebsch)>1.e-8 ){
	    accept=true;
            if( (lampbar-lamp) == Spin(0) ) spinProjection0=true;
            else if( (lampbar-lamp) == Spin(1) ) spinProjection1=true;	    
	  }
	}
      }     
    }
    if(accept){
      _pbarpLSs.push_back(LScandidates.at(ls));
      if( S==1 ){
	_isTriplet=true;
	if(spinProjection0) _isTriplet0=true;
	if(spinProjection1) _isTriplet1=true;
      }
      else if( S==0 ) _isSinglet=true;
    }        
  }
}

void pbarpState::print(std::ostream& os) const{
  os << "JPC: " << _jpc->name();
  os << "*** allowed LS combinations for the pbar p system *** " << std::endl;
  std::vector<LS>::const_iterator itLS;
 for ( itLS=_pbarpLSs.begin(); itLS!=_pbarpLSs.end(); ++itLS){
   os << "L: " << itLS->L << "\tS: " << itLS->L << std::endl; 
 }
 if(_isSinglet) os << "is singlet state with spin projection 0" << std::endl;
 if(_isTriplet){ 
   if(_isTriplet0 && !_isTriplet1) os << " triplet state with spin projection 0" << std::endl;
   if(!_isTriplet0 && _isTriplet1) os << " triplet state with spin projections +-1" << std::endl;
   if(_isTriplet0 && _isTriplet1)  os << " triplet state with spin projections 0 and +-1" << std::endl;;   
 }  
}

