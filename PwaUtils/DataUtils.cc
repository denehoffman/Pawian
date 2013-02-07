// DataUtils class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include "PwaUtils/DataUtils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"



void validJPCLS(boost::shared_ptr<const jpcRes> motherRes, boost::shared_ptr<const jpcRes> daughterRes1, boost::shared_ptr<const jpcRes> daughterRes2, std::vector< boost::shared_ptr<const JPCLS> >& theJPCLSVec)
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
    Spin L= LSs[ls].L; 
    Spin S= LSs[ls].S;
    boost::shared_ptr<const JPCLS> tmpJPCLS(new JPCLS(motherRes, L, S));
    theJPCLSVec.push_back(tmpJPCLS);
  }
}

void validJPCLS(boost::shared_ptr<const jpcRes> motherRes, Particle* daughter1, Particle* daughter2, std::vector< boost::shared_ptr<const JPCLS> >& theJPCLSVec){

  // first: check C-parity
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

  bool identicalDaughters=false;
  if( (*daughter1)==(*daughter2)) identicalDaughters=true;
  vector<LS> LSs=GetValidLS(motherRes->J, motherRes->P, daughter1->J(), daughter1->theParity(), daughter2->J(), daughter2->theParity());

  int num_LS = (int) LSs.size();
  for(int ls = 0; ls < num_LS; ls++){
    Spin L= LSs[ls].L; 
    Spin S= LSs[ls].S;
    bool LplusSeven=false;
    int LplusS=L+S;
    if( LplusS% 2 == 0) LplusSeven=true;
    if(!identicalDaughters || (identicalDaughters && LplusSeven)){ 
      boost::shared_ptr<const JPCLS> tmpJPCLS(new JPCLS(motherRes, L, S));
      theJPCLSVec.push_back(tmpJPCLS);
    }
  }
}


void validJPClamlam(boost::shared_ptr<const jpcRes> motherRes, Particle* daughter1, Particle* daughter2, std::vector< boost::shared_ptr<const JPClamlam> >& theJPClamlamVec){

  std::vector< boost::shared_ptr<const JPCLS> > currentJPCLSDecAmps;
  boost::shared_ptr<const jpcRes> daughterRes1=getJPCPtr(daughter1);
  boost::shared_ptr<const jpcRes> daughterRes2=getJPCPtr(daughter2);

  validJPCLS(motherRes, daughterRes1, daughterRes2, currentJPCLSDecAmps);
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

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator it;
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
	std::vector< boost::shared_ptr<const JPClamlam> >::iterator it;
	for (it=theJPClamlamVec.begin(); it!=theJPClamlamVec.end(); ++it){
	  if( ((*it)->lam1==lam2 && (*it)->lam2==lam1) ||  ((*it)->lam1==-lam2 && (*it)->lam2==-lam1) ){
	    fillIt=false;
	  }
	} 
      }

      if(fillIt){
	boost::shared_ptr<const JPClamlam> tmpJPClamlam(new JPClamlam(motherRes, lam1, lam2, 1.));
	theJPClamlamVec.push_back(tmpJPClamlam);
      }
    }
  }  
}


boost::shared_ptr<jpcRes> getJPCPtr(Particle* theParticle){

  boost::shared_ptr<jpcRes> result(new jpcRes((int) theParticle->J(), theParticle->theParity(), theParticle->theCParity()));

  return result;  
}

