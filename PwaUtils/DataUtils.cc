#include "PwaUtils/DataUtils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Utils/PawianCollectionUtils.hh"
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


boost::shared_ptr<jpcRes> getJPCPtr(Particle* theParticle){

  boost::shared_ptr<jpcRes> result(new jpcRes((int) theParticle->J(), theParticle->theParity(), theParticle->theCParity()));

  return result;  
}
