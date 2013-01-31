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
  // std::cout << "daughter1->theParity()\t" << daughter1->theParity() << "\tdaughter2->theParity()\t" << daughter1->theParity() <<std::endl; 
  // std::cout << "Jmother=\t" << Jmother << "\tJdaughter1=\t" << Jdaughter1 << "\tJdaughter2=\t" << Jdaughter2 <<"\tparityFactor=\t" << parityFactor << "\n" <<std::endl; 
  bool isDaughter1Photon=false;
  if(daughter1->name()=="photon") isDaughter1Photon=true;
  bool isDaughter2Photon=false;
  if(daughter2->name()=="photon") isDaughter2Photon=true;

  for (Spin lam1=-Jdaughter1; lam1<=Jdaughter1; ++lam1){
    if(isDaughter1Photon && lam1==0) continue;
    for (Spin lam2=lam1; lam2<=Jdaughter2; ++lam2){
      if(lam1==0 && lam2==0 && parityFactor <0) continue; //parity conservation 
      if (lam2<-Jdaughter2) continue; 
      if(isDaughter2Photon && lam2==0) continue;
      Spin lambda=lam1-lam2;
      if (fabs(lambda)>Smax) continue;
      boost::shared_ptr<const JPClamlam> tmpJPClamlam(new JPClamlam(motherRes, lam1, lam2));
      theJPClamlamVec.push_back(tmpJPClamlam);
    }
  }  
}


boost::shared_ptr<jpcRes> getJPCPtr(Particle* theParticle){

  boost::shared_ptr<jpcRes> result(new jpcRes((int) theParticle->J(), theParticle->theParity(), theParticle->theCParity()));

  return result;  
}

