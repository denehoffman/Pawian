#include "Examples/MATpbarpToOmegaPi/AbsStates.hh"
#include "Utils/MathUtils.hh"

AbsStates::AbsStates(){
}


AbsStates::~AbsStates(){
}

std::vector< boost::shared_ptr<const JPCLSM> > AbsStates::extractJPCLSMStates(std::vector< boost::shared_ptr<const jpcRes> >& theJPCStates) const{

  std::vector< boost::shared_ptr<const JPCLSM> > result;

  std::vector< boost::shared_ptr<const jpcRes> >::const_iterator itJPC;
  std::vector< boost::shared_ptr<const JPCLSM> >::const_iterator itAllStates;

  for ( itJPC=theJPCStates.begin(); itJPC!=theJPCStates.end(); ++itJPC){
    const jpcRes* jpcRequest=(*itJPC).get();

    for ( itAllStates=_allStates.begin(); itAllStates!=_allStates.end(); ++itAllStates){
      const JPCLSM* jpcsmCurrent=(*itAllStates).get(); 
      if (( *jpcsmCurrent) ==  (*jpcRequest)){
	result.push_back(*itAllStates);
	continue;
      }
    }
  }
  return result;
}

std::vector< boost::shared_ptr<const JPCSM> > AbsStates::extractJPCSMStates(std::vector< boost::shared_ptr<const jpcRes> >& theJPCStates) const{

  std::vector< boost::shared_ptr<const JPCSM> > result;

  std::vector< boost::shared_ptr<const jpcRes> >::const_iterator itJPC;
  std::vector< boost::shared_ptr<const JPCSM> >::const_iterator itAllStates;

  for ( itJPC=theJPCStates.begin(); itJPC!=theJPCStates.end(); ++itJPC){
    const jpcRes* jpcRequest=(*itJPC).get();

    for ( itAllStates=_allJPCSM.begin(); itAllStates!=_allJPCSM.end(); ++itAllStates){
      const JPCSM* jpcsmCurrent=(*itAllStates).get(); 
      if (( *jpcsmCurrent) ==  (*jpcRequest)){
	result.push_back(*itAllStates);
	continue;
      }
    }
  }
  return result;

}

std::vector< boost::shared_ptr<const JPCLS> > AbsStates::extractJPCLSStates(boost::shared_ptr<const jpcRes> theJPCState) const{
  std::vector< boost::shared_ptr<const JPCLS> > result;

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator it;
  for ( it=_allJPCLS.begin(); it!=_allJPCLS.end(); ++it){
//     const jpcRes* currentJPC= (jpcRes*) (*it).get(); 
//     const jpcRes* currentJPC=dynamic_cast<const jpcRes*> ((*it).get());

//     if ( *(currentJPC) == *(theJPCState.get())) result.push_back( (*it) );
    if (*((*it).get())==*(theJPCState.get())) result.push_back( (*it) );
  }
 
  return result;
}

void AbsStates::print(std::ostream& os) const{

  std::vector< boost::shared_ptr<const JPCLSM> >::const_iterator it;
   for ( it=_allStates.begin(); it!=_allStates.end(); ++it){
    if (0!= (*it)){
      (*it)->print(os); 
// 	os <<"\tL=" << (*it)->lsm.L <<"\tS=" << (*it)->lsm.S <<"\tlambda=" << (*it)->lsm.M
// 	   <<"\tClebschGordan=" << (*it)->ClebschG 
// 	   << std::endl;
      
    }
   }
}

