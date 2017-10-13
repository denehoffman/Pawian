//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// gamgamStates class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#include "gamgamUtils/gamgamStates.hh"
#include "Utils/MathUtils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"

gamgamStates::gamgamStates():
  AbsStates(),
  _jmax(10),
  _gammaJPC(1,-1,-1)
{
  //  GlobalEnv::instance()->setup();
  _gammaParticle=GlobalEnv::instance()->particleTable()->particle("photon");
   calcStates();
}


gamgamStates::gamgamStates(int jmax):
  AbsStates(),
  _jmax(jmax),
  _gammaJPC(1,-1,-1)
{
  //   GlobalEnv::instance()->setup();
   _gammaParticle=GlobalEnv::instance()->particleTable()->particle("photon");
    calcStates();
}

gamgamStates::~gamgamStates(){
}

std::vector< std::shared_ptr<const JPCLS> > gamgamStates::jpcLSStates(Spin lamGam1, Spin lamGam2) const{
  if(lamGam1==-1 && lamGam2==-1) return _JPCLSLam1m1Lam2m1_States;
  else if(lamGam1==-1 && lamGam2==1) return _JPCLSLam1m1Lam2p1_States;
  else if(lamGam1==1 && lamGam2==-1) return _JPCLSLam1p1Lam2m1_States;
  else if(lamGam1==1 && lamGam2==1) return _JPCLSLam1p1Lam2p1_States;
  else{
    Alert << "lamGam1= " << lamGam1 << " and lamGam2= " << lamGam2 << " are not allowed for two initial real photons!!!" << endmsg;
    exit(1);
  } 
}

std::vector< std::shared_ptr<const JPClamlam> > gamgamStates::jpcLamLamStates(Spin lamGam1, Spin lamGam2) const{
  if(lamGam1==-1 && lamGam2==-1) return _JPClamlamLam1m1Lam2m1_States;
  else if(lamGam1==-1 && lamGam2==1) return _JPClamlamLam1m1Lam2p1_States;
  else if(lamGam1==1 && lamGam2==-1) return _JPClamlamLam1p1Lam2m1_States;
  else if(lamGam1==1 && lamGam2==1) return _JPClamlamLam1p1Lam2p1_States;
  else{
    Alert << "lamGam1= " << lamGam1 << " and lamGam2= " << lamGam2 << " are not allowed for two initial real photons!!!" << endmsg;
    exit(1);
  } 
}

bool gamgamStates::calcStates(){
  std::cout << "jmax: " << _jmax << std::endl;
  for (int j=0; j<=_jmax; j++){
    for(int parity=-1; parity<=1; parity+=2){
    int cparity(1);
    std::shared_ptr<const jpcRes> jpcPtr(new jpcRes(j, parity, cparity));
    for (int iso=0; iso<=2; iso++){ //isospin
      for (int gp=-1; gp<=1; gp+=2){ //gparity
	//first check g-parity
	if (gp != _gammaParticle->theGParity()*_gammaParticle->theGParity() ) continue;
	std::shared_ptr<const IGJPC> currentIGJPCPtr(new IGJPC(jpcPtr, iso, gp, 0.));
	std::vector< std::shared_ptr<const LScomb> > theLSVec;
	bool useCParity(true);
	bool useIsospin(true);
	validLS( jpcPtr, _gammaParticle, _gammaParticle, theLSVec, useCParity, gp, useIsospin);
	std::vector< std::shared_ptr<const LScomb> >::iterator it;
	for(it=theLSVec.begin(); it!=theLSVec.end();++it){
	  //look for JPClamlam states
	  for(Spin lam1=-1; lam1<=1; lam1+=2){
	    for(Spin lam2=-1; lam2<=1; lam2+=2){
	      double preFactorLS=sqrt(2.*(*it)->L+1)*Clebsch((*it)->S, lam1-lam2, (*it)->L, 0, jpcPtr->J,  lam1-lam2);
	      if( fabs(preFactorLS) > 1.e-10){
		std::shared_ptr<const JPCLS> currentJPCLSPtr(new JPCLS(jpcPtr, (*it)->L, (*it)->S, preFactorLS)); 
		double cgFactorLamLam = Clebsch(1, lam1, 1, -lam2, (*it)->S,  lam1-lam2);
		double parityFactorLamLam = jpcPtr->P*pow(-1., jpcPtr->J);
		if(fabs(cgFactorLamLam) > 1.e-10){
		  fillVec(currentJPCLSPtr,_JPCLS_States);
		  fillVec(currentIGJPCPtr, _allIGjpcRes);
		  std::shared_ptr<const JPClamlam> currentJPClamlamPtr(new JPClamlam(jpcPtr, lam1, lam2, parityFactorLamLam));
		  fillVec(currentJPClamlamPtr, _JPClamlam_States);
		  fillVec(jpcPtr, _alljpcRes);
		  if(lam1==-1 && lam2==-1){
		    fillVec(currentJPCLSPtr, _JPCLSLam1m1Lam2m1_States);
		    fillVec(currentJPClamlamPtr, _JPClamlamLam1m1Lam2m1_States);
		  }
		  else if(lam1==-1 && lam2==1){
		    fillVec(currentJPCLSPtr, _JPCLSLam1m1Lam2p1_States);
		    fillVec(currentJPClamlamPtr, _JPClamlamLam1m1Lam2p1_States);
		  }
		  else if(lam1==1 && lam2==-1){
		    fillVec(currentJPCLSPtr, _JPCLSLam1p1Lam2m1_States);
		    fillVec(currentJPClamlamPtr, _JPClamlamLam1p1Lam2m1_States);
		  }
		  else if(lam1==1 && lam2==1){
		    fillVec(currentJPCLSPtr, _JPCLSLam1p1Lam2p1_States);
		    fillVec(currentJPClamlamPtr, _JPClamlamLam1p1Lam2p1_States);
		  } 
		}
	      }
	    }
	  }
	}
      }
    }
    }     
    
  }
  
  return true;
}


void gamgamStates::print(std::ostream& os) const{
   std::vector< std::shared_ptr<const JPCLS > >::const_iterator itJPCLS;
   std::vector< std::shared_ptr<const jpcRes > >::const_iterator itjpcRes;
   std::vector< std::shared_ptr<const IGJPC > >::const_iterator itIGJPC;
   std::vector< std::shared_ptr<const JPClamlam > >::const_iterator itJPCLamLam;
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

 os << "\n******** JPCLS states  ************ " << std::endl;
 for ( itJPCLS=_JPCLS_States.begin(); itJPCLS!=_JPCLS_States.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << std::endl;
  }

 os << "\n******** JPCLamLam states  ************ " << std::endl;
 for ( itJPCLamLam=_JPClamlam_States.begin(); itJPCLamLam!=_JPClamlam_States.end(); ++itJPCLamLam){
   (*itJPCLamLam)->print(os);
   os << std::endl;
 }

 for (Spin lam1=-1; lam1<=1; lam1+=2){
   for (Spin lam2=-1; lam2<=1; lam2+=2){
     os << "\n******** JPCLS states  for lamGam1= " << lam1 << " and lamGam2= " << lam2 <<" ===> lamX= " << lam1-lam2 << " ************ " << std::endl;
     std::vector< std::shared_ptr<const JPCLS> > currentJPCLSStates=jpcLSStates(lam1, lam2);
     for ( itJPCLS=currentJPCLSStates.begin(); itJPCLS!=currentJPCLSStates.end(); ++itJPCLS){
       (*itJPCLS)->print(os);
       os << std::endl;
     }
   }
 }

 for (Spin lam1=-1; lam1<=1; lam1+=2){
   for (Spin lam2=-1; lam2<=1; lam2+=2){
     os << "\n******** JPClamlam states  for lamGam1= " << lam1 << " and lamGam2= " << lam2 <<" ===> lamX= " << lam1-lam2 << " ************ " << std::endl;
     std::vector< std::shared_ptr<const JPClamlam> > currentJPClamalamStates=jpcLamLamStates(lam1, lam2);
     for ( itJPCLamLam=currentJPClamalamStates.begin(); itJPCLamLam!=currentJPClamalamStates.end(); ++itJPCLamLam){
       (*itJPCLamLam)->print(os);
       os << std::endl;
     }
   }
 }
}

