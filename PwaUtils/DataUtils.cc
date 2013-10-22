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

// DataUtils class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include "PwaUtils/DataUtils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"



void validJPCLS(std::shared_ptr<const jpcRes> motherRes, std::shared_ptr<const jpcRes> daughterRes1, std::shared_ptr<const jpcRes> daughterRes2, std::vector< std::shared_ptr<const JPCLS> >& theJPCLSVec)
{
  std::vector< std::shared_ptr<const LScomb> > LSVec;
  validLS(motherRes, daughterRes1, daughterRes2, LSVec);

  std::vector< std::shared_ptr<const LScomb> >::iterator it;

  for(it=LSVec.begin(); it!=LSVec.end(); ++it){
    std::shared_ptr<const JPCLS> tmpJPCLS(new JPCLS(motherRes, (*it)->L, (*it)->S));
    theJPCLSVec.push_back(tmpJPCLS);
  }
}

void validLS(std::shared_ptr<const jpcRes> motherRes, std::shared_ptr<const jpcRes> daughterRes1, std::shared_ptr<const jpcRes> daughterRes2, std::vector< std::shared_ptr<const LScomb> >& theLSVec)
{
  // first: check C-parity
  if ( motherRes->C != daughterRes1->C*daughterRes2->C){
    Warning << "C-Parity not valid for the reaction: JPC= " 
            << motherRes->J << " " << motherRes->P << " " << motherRes->C
            << " --> "
            << " JPC= " << daughterRes1->J << " " << daughterRes1->P << " " << daughterRes1->C
            << " and "
            << " JPC= " << daughterRes2->J << " " << daughterRes2->P << " " << daughterRes2->C
            ;  // << endmsg;
    if( fabs(motherRes->C)==1 && fabs(daughterRes1->C)==1 && fabs(daughterRes2->C)==1) return; 
  }

 
  vector<LS> LSs=GetValidLS(motherRes->J, motherRes->P, daughterRes1->J, daughterRes1->P, daughterRes2->J, daughterRes2->P);

  int num_LS = (int) LSs.size();

  for(int ls = 0; ls < num_LS; ls++){
    int L= LSs[ls].L; 
    Spin S= LSs[ls].S;
    std::shared_ptr<const LScomb> tmpLS(new LScomb(L, S));
    theLSVec.push_back(tmpLS);
  }

  if(LSs.size()==0) Info << "size for decay " << motherRes->name() << " to " << daughterRes1->name() << " and " << daughterRes2->name()
                         << " =0!!!" << endmsg;
}


void validJPCLS(std::shared_ptr<const jpcRes> motherRes, Particle* daughter1, Particle* daughter2, std::vector< std::shared_ptr<const JPCLS> >& theJPCLSVec, bool useCParity, int gParityMother, bool useIsospin){

  std::vector< std::shared_ptr<const LScomb> > LSVec;
  validLS(motherRes, daughter1, daughter2, LSVec, useCParity, gParityMother, useIsospin);

  std::vector< std::shared_ptr<const LScomb> >::iterator it;

  for(it=LSVec.begin(); it!=LSVec.end(); ++it){
    std::shared_ptr<const JPCLS> tmpJPCLS(new JPCLS(motherRes, (*it)->L, (*it)->S));
    theJPCLSVec.push_back(tmpJPCLS);
  }
}

void validLS(std::shared_ptr<const jpcRes> motherRes, Particle* daughter1, Particle* daughter2, std::vector< std::shared_ptr<const LScomb> >& theLSVec, bool useCParity, int gParityMother, bool useIsospin){
  // first: check C-parity
  if (useCParity){
    if ( motherRes->C != daughter1->theCParity()*daughter2->theCParity()){
      Warning << "C-Parity not valid for the reaction: JPC= " 
              << motherRes->J << " " << motherRes->P << " " << motherRes->C
              << " --> "
              << " JPC= " << daughter1->J() << " " << daughter1->theParity() << " " << daughter1->theCParity()
              << " and "
              << " JPC= " << daughter2->J() << " " << daughter2->theParity() << " " << daughter2->theCParity()
        ;  // << endmsg;
      if( fabs(motherRes->C)==1 && fabs(daughter1->theCParity())==1 && fabs(daughter2->theCParity())==1) return; 
    }
  }

  if(useIsospin){
    // second: check G-parity
    if (gParityMother != daughter1->theGParity()*daughter2->theGParity() ){
      Warning << "G-Parity not valid for:" 
              << gParityMother 
              << " --> "
              << daughter1->theGParity() << " * " << daughter2->theGParity() << " " 
              << endmsg;
      
      if( fabs(gParityMother)==1 && fabs(daughter1->theGParity())==1 && fabs(daughter2->theGParity())==1) return; 
    }
  }

  bool identicalDaughters=false;
  if( (*daughter1)==(*daughter2)) identicalDaughters=true;
  vector<LS> LSs=GetValidLS(motherRes->J, motherRes->P, daughter1->J(), daughter1->theParity(), daughter2->J(), daughter2->theParity());

  int num_LS = (int) LSs.size();
  for(int ls = 0; ls < num_LS; ls++){
    int L= LSs[ls].L; 
    Spin S= LSs[ls].S;
    bool LplusSeven=false;
    int LplusS=L+S;
    if( LplusS% 2 == 0) LplusSeven=true;
    if(!identicalDaughters || (identicalDaughters && LplusSeven)){ 
      std::shared_ptr<const LScomb> tmpLS(new LScomb(L, S));
      theLSVec.push_back(tmpLS);
    }
  }
}


void validJPClamlam(std::shared_ptr<const jpcRes> motherRes, Particle* daughter1, Particle* daughter2, std::vector< std::shared_ptr<const JPClamlam> >& theJPClamlamVec, bool useCParity, int gParityMother, bool useIsospin){
  std::vector< std::shared_ptr<const JPCLS> > currentJPCLSDecAmps;
  std::shared_ptr<const jpcRes> daughterRes1=getJPCPtr(daughter1);
  std::shared_ptr<const jpcRes> daughterRes2=getJPCPtr(daughter2);

  //  validJPCLS(motherRes, daughterRes1, daughterRes2, currentJPCLSDecAmps, useCParity, gParityMother, useIsospin);
  validJPCLS(motherRes, daughter1, daughter2, currentJPCLSDecAmps, useCParity, gParityMother, useIsospin);
  if(currentJPCLSDecAmps.size()==0){
    Warning << "decay JPC= " 
	  << motherRes->J << " " << motherRes->P << " " << motherRes->C
	  << " --> "
	  << " JPC= " << daughterRes1->J << " " << daughterRes1->P << " " << daughterRes1->C
	  << " and "
	  << " JPC= " << daughterRes2->J << " " << daughterRes2->P << " " << daughterRes2->C 
	  << " not allowed!!!" 
	  << endmsg;
    return;
  }

  bool identicalDaughters=false;
  if( (*daughter1)==(*daughter2)) identicalDaughters=true;

  Spin Smax=0;

  std::vector< std::shared_ptr<const JPCLS> >::const_iterator it;
  std::vector<Spin> itSpin;
  for (it=currentJPCLSDecAmps.begin(); it!=currentJPCLSDecAmps.end(); ++it){
    Spin currentSpin=(*it)->S;
    if (currentSpin>Smax) Smax=currentSpin;
  }

  Spin Jmother=motherRes->J;
  Spin Jdaughter1=daughterRes1->J;
  Spin Jdaughter2=daughterRes2->J;

  double parityFactor = motherRes->P*daughter1->theParity()*daughter2->theParity()*pow(-1., Jmother - Jdaughter1 - Jdaughter2 );



  bool isDaughter1Photon=false;
  if(daughter1->name()=="photon") isDaughter1Photon=true;
  bool isDaughter2Photon=false;
  if(daughter2->name()=="photon") isDaughter2Photon=true;

  for (Spin lam1=Jdaughter1; lam1>=-Jdaughter1; --lam1){
    if(isDaughter1Photon && fabs(lam1)!=1) continue;
    for (Spin lam2=lam1; lam2>=-Jdaughter2; --lam2){
      if (fabs(lam2)>Jdaughter2) continue;
      if(lam1==0 && lam2==0 && parityFactor <0) continue; //parity conservation
      if( identicalDaughters && parityFactor <0 && lam1==-lam2 && lam2==-lam1) continue; //parity conservation + identical particles 
      if(isDaughter2Photon && fabs(lam2)!=1) continue;
      if(lam1==lam2 && lam1<0) continue;
      Spin lambda=lam1-lam2;
      if (fabs(lambda)>Smax || fabs(lambda)>Jmother) continue;
      bool fillIt=true;
      if( identicalDaughters ){
	//check if lam2 lam1 combination already exists
	std::vector< std::shared_ptr<const JPClamlam> >::iterator it;
	for (it=theJPClamlamVec.begin(); it!=theJPClamlamVec.end(); ++it){
	  if( ((*it)->lam1==lam2 && (*it)->lam2==lam1) ||  ((*it)->lam1==-lam2 && (*it)->lam2==-lam1) ){
	    fillIt=false;
	  }
	} 
      }

      if(fillIt){
	std::shared_ptr<const JPClamlam> tmpJPClamlam(new JPClamlam(motherRes, lam1, lam2, 1.));
	theJPClamlamVec.push_back(tmpJPClamlam);
      }
    }
  }  
}


std::shared_ptr<jpcRes> getJPCPtr(Particle* theParticle){

  std::shared_ptr<jpcRes> result(new jpcRes((int) theParticle->J(), theParticle->theParity(), theParticle->theCParity()));

  return result;  
}

std::shared_ptr<IGJPC> getIGJPCPtr(Particle* theParticle){
  Spin isoSpin=Spin(theParticle->twoIso3(), 2);
  std::shared_ptr<IGJPC> result(new IGJPC((int) theParticle->J(), theParticle->theParity(), theParticle->theCParity(), isoSpin, theParticle->theGParity()));

  return result;  
}

