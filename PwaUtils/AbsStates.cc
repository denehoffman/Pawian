#include "PwaUtils/AbsStates.hh"

#include "Utils/MathUtils.hh"
#include "ErrLogger/ErrLineLog.hh"

 

AbsStates::AbsStates(){
}


AbsStates::~AbsStates(){
}

std::vector< boost::shared_ptr<const JPCSML> > AbsStates::extractJPCSMLStates(std::vector< boost::shared_ptr<const jpcRes> >& theJPCStates){

  std::vector< boost::shared_ptr<const JPCSML> > result;

  std::vector< boost::shared_ptr<const jpcRes> >::const_iterator itJPC;
  std::vector< boost::shared_ptr<const JPCSML> >::const_iterator itAllStates;

  for ( itJPC=theJPCStates.begin(); itJPC!=theJPCStates.end(); ++itJPC){
    const jpcRes* jpcRequest=(*itJPC).get();

    for ( itAllStates=_allStates.begin(); itAllStates!=_allStates.end(); ++itAllStates){
      const jpcRes* jpcCurrent=(*itAllStates)->jpc.get(); 
      if (( *jpcCurrent) ==  (*jpcRequest)){
	result.push_back(*itAllStates);
	continue;
      }
    }
  }
  return result;
}

std::vector< boost::shared_ptr<const JPCSM> > AbsStates::extractJPCSMStates(std::vector< boost::shared_ptr<const jpcRes> >& theJPCStates){

  std::vector< boost::shared_ptr<const JPCSM> > result;

  std::vector< boost::shared_ptr<const jpcRes> >::const_iterator itJPC;
  std::vector< boost::shared_ptr<const JPCSM> >::const_iterator itAllStates;

  for ( itJPC=theJPCStates.begin(); itJPC!=theJPCStates.end(); ++itJPC){
    const jpcRes* jpcRequest=(*itJPC).get();

    for ( itAllStates=_allJPCSM.begin(); itAllStates!=_allJPCSM.end(); ++itAllStates){
      const jpcRes* jpcCurrent=(*itAllStates)->jpc.get(); 
      if (( *jpcCurrent) ==  (*jpcRequest)){
	result.push_back(*itAllStates);
	continue;
      }
    }
  }
  return result;

}



void AbsStates::print(std::ostream& os) const{

  std::vector< boost::shared_ptr<const JPCSML> >::const_iterator it;
   for ( it=_allStates.begin(); it!=_allStates.end(); ++it){
    if (0!= (*it)){
      (*it)->jpc->print(os); 
	os <<"\tL=" << (*it)->lsm.L <<"\tS=" << (*it)->lsm.S <<"\tlambda=" << (*it)->lsm.M
	   <<"\tClebschGordan=" << (*it)->ClebschG 
	   << std::endl;
      
    }
   }
}

