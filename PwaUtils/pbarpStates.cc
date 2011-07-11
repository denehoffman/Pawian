#include "PwaUtils/pbarpStates.hh"
#include "Utils/MathUtils.hh"
#include "ErrLogger/ErrLogger.hh"

pbarpStates::pbarpStates():
  AbsStates(),
  _jmax(10),
  _pbarJPC(0.5, -1),
  _pJPC(0.5, 1)
{
  calcStates();
}


pbarpStates::pbarpStates(int jmax):
  AbsStates(),
  _jmax(jmax),
  _pbarJPC(0.5, -1),
  _pJPC(0.5, 1)
{
  calcStates();
}

pbarpStates::~pbarpStates(){
}

bool pbarpStates::calcStates(){
  for (int j=0; j<=_jmax; j++)
    {
      for (int p=-1; p<=1; p+=2){
       vector<LS> motherLSs=GetValidLS(j, p, _pbarJPC.J,_pbarJPC.P, _pJPC.J, _pJPC.P);
       int num_ls = (int) motherLSs.size();
       DebugMsg << "valid LS combinations for JP(pbar p)=" << j <<" " << p << endmsg;
       
       for(int ls = 0; ls < num_ls; ls++){

	 Spin L= motherLSs[ls].L;
	 Spin S= motherLSs[ls].S;
	 int cparity(pow(-1,int(L+S)));
         DebugMsg << "L=" << L <<" S=" << S << " ==> C=" << cparity << endmsg;
	 boost::shared_ptr<const jpcRes> jpcPtr(new jpcRes(j,p,cparity));
 
	 std::vector<Spin> theMs;
        
	 for(Spin M = -S; M <= S; M++){
           double Clebschg=Clebsch(L,0,S,M, j,M);
	   DebugMsg << "Clebsch(L,0,S,M=" << M << ", j,M=" << M << "): " << Clebschg << endmsg;
           if (fabs(Clebschg)>1e-8) theMs.push_back(M);
	 }


	 if (theMs.size()>0){
	   fillVec(jpcPtr, _jpcStates);	   

	   boost::shared_ptr<const JPCLS> tmpJPCLS(new JPCLS(jpcPtr, L, S) );
	   _allJPCLS.push_back(tmpJPCLS);
	   
           std::vector<Spin>::const_iterator itM;
	   for ( itM=theMs.begin(); itM!=theMs.end(); ++itM){
	     boost::shared_ptr<const JPCSM> tmpJPCSM( new JPCSM(jpcPtr, S, (*itM) ));
	     _allJPCSM .push_back(tmpJPCSM);
	     
	     double Clebschg1=Clebsch(L,0,S,(*itM), j, (*itM));
	     boost::shared_ptr<const JPCLSM> tmpJPCLSM(new JPCLSM(tmpJPCLS, (*itM), Clebschg1) );
	     _allStates.push_back(tmpJPCLSM);
             if (S==0) fillVec(jpcPtr, _singletStates);
	     else if (S==1){
	       if ((*itM)==0) fillVec(jpcPtr, _tripletM0States);
               else if ((*itM)==1) fillVec(jpcPtr, _tripletMp1States);
               else if ((*itM)==-1) fillVec(jpcPtr, _tripletMm1States);
	       else {
		 Alert << "pbar p state with S=" << S << " and M="<< (*itM) <<" cannot exitst!!!" 
		       << cparity << endmsg;
		 exit(1);
	       }
	     }
	     else {
	       Alert << "pbar p state with S=" << S << " cannot exitst!!!" << cparity << endmsg;
	       exit(1);
	     }
	   }             
	 }
        
       }
      }
       DebugMsg << "\n" << endmsg;
    }

  return true;
}

void pbarpStates::fillVec(boost::shared_ptr<const jpcRes> currentRes, 
			  std::vector< boost::shared_ptr<const jpcRes> >& theVec){
  const jpcRes* jpcCurrent=currentRes.get();
  if (0==jpcCurrent) {
    Alert << "shared object containss 0 pointer!!!" << endmsg;
    exit(1);
  }
  std::vector< boost::shared_ptr<const jpcRes> >::const_iterator it;
  bool found=false;
  for (it=theVec.begin(); it!=theVec.end(); ++it){
    const jpcRes* jpcIt=(*it).get();
    if (0==jpcIt) {
      Alert << "shared object containss 0 pointer!!!" << endmsg;
      exit(1);
    }
    if ( (*jpcIt) == (*jpcCurrent) ){
      found=true;
      continue;
    }
  }
  if (!found) theVec.push_back(currentRes);  

}

std::vector< boost::shared_ptr<const jpcRes> > pbarpStates::extractSingletStates(std::vector< boost::shared_ptr<const jpcRes> >& theJPCStates) const{
  return extractJPCStates(theJPCStates,_singletStates);
}

std::vector< boost::shared_ptr<const jpcRes> > pbarpStates::extractTripletM0States(std::vector< boost::shared_ptr<const jpcRes> >& theJPCStates) const{
   return extractJPCStates(theJPCStates, _tripletM0States);
}

std::vector< boost::shared_ptr<const jpcRes> > pbarpStates::extractTripletM1States(std::vector< boost::shared_ptr<const jpcRes> >& theJPCStates) const{
   return extractJPCStates(theJPCStates, _tripletMp1States);
}

std::vector< boost::shared_ptr<const jpcRes> > pbarpStates::extractJPCStates(std::vector< boost::shared_ptr<const jpcRes> > theJPCRequests, std::vector< boost::shared_ptr<const jpcRes> > jpcStatesAll) const{

  std::vector< boost::shared_ptr<const jpcRes> > result;

  std::vector< boost::shared_ptr<const jpcRes> >::const_iterator itJPC;
  std::vector< boost::shared_ptr<const jpcRes> >::const_iterator itAllStates;
  for ( itJPC=theJPCRequests.begin(); itJPC!=theJPCRequests.end(); ++itJPC){
    const jpcRes* jpcRequest=(*itJPC).get();

    for ( itAllStates=jpcStatesAll.begin(); itAllStates!=jpcStatesAll.end(); ++itAllStates){
      const jpcRes* jpcCurrent=(*itAllStates).get();
      if (( *jpcCurrent) ==  (*jpcRequest)){
	result.push_back(*itAllStates);
	continue;
      } 
    }
  }
  return result;

}

std::vector< boost::shared_ptr<const JPCLSM> > pbarpStates::extractJPCLSMStates(std::vector< boost::shared_ptr<const jpcRes> >& theJPCStates) const{

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

std::vector< boost::shared_ptr<const JPCSM> > pbarpStates::extractJPCSMStates(std::vector< boost::shared_ptr<const jpcRes> >& theJPCStates) const{

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

std::vector< boost::shared_ptr<const JPCLS> > pbarpStates::extractJPCLSStates(boost::shared_ptr<const jpcRes> theJPCState) const{
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

void pbarpStates::print(std::ostream& os) const{
  os << "initital states of the pbar p annihilation for Jmax = " << _jmax << " are: " << std::endl; 
  AbsStates::print(os);

}

