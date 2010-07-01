#include "PwaUtils/pbarpStates.hh"

#include "Utils/MathUtils.hh"

#include "ErrLogger/ErrLineLog.hh"


 

pbarpStates::pbarpStates():
  _jmax(10),
  _pbarJPC(0.5, -1),
  _pJPC(0.5, 1),
  AbsStates()
{
  calcJPCs();
}


pbarpStates::pbarpStates(int jmax):
  _jmax(jmax),
  _pbarJPC(0.5, -1),
  _pJPC(0.5, 1),
  AbsStates()
{
  calcJPCs();
}

pbarpStates::~pbarpStates(){
}

bool pbarpStates::calcJPCs(){
  for (int j=0; j<=_jmax; j++)
    {
      for (int p=-1; p<=1; p+=2){
       vector<LS> motherLSs=GetValidLS(j, p, _pbarJPC.J,_pbarJPC.P, _pJPC.J, _pJPC.P);
       int num_ls = (int) motherLSs.size();
       ErrMsg(debugging) << "valid LS combinations for JP(pbar p)=" << j <<" " << p << endmsg;
       
       for(int ls = 0; ls < num_ls; ls++){

	 Spin L= motherLSs[ls].L;
	 Spin S= motherLSs[ls].S;
	 int cparity(pow(-1,int(L+S)));
         ErrMsg(debugging) << "L=" << L <<" S=" << S << " ==> C=" << cparity << endmsg;
	 boost::shared_ptr<const jpcRes> jpcPtr(new jpcRes(j,p,cparity));
 
	 std::vector<Spin> theMs;
        
	 for(Spin M = -S; M <= S; M++){
           double Clebschg=Clebsch(L,0,S,M, j,M);
	   ErrMsg(debugging) << "Clebsch(L,0,S,M=" << M << ", j,M=" << M << "): " << Clebschg << endmsg;
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
	       else ErrMsg(fatal) << "pbar p state with S=" << S << " and M="<< (*itM) <<" cannot exitst!!!" << cparity << endmsg;
	     }
	     else ErrMsg(fatal) << "pbar p state with S=" << S << " cannot exitst!!!" << cparity << endmsg;
	   }             
	 }
        
       }
      }
       ErrMsg(debugging) << "\n" << endmsg;
    }
}

void pbarpStates::fillVec(boost::shared_ptr<const jpcRes> currentRes, std::vector< boost::shared_ptr<const jpcRes> >& theVec){
  const jpcRes* jpcCurrent=currentRes.get();
  if (0==jpcCurrent) ErrMsg(fatal) << "shared object containss 0 pointer!!!" << endmsg;
  std::vector< boost::shared_ptr<const jpcRes> >::const_iterator it;
  bool found=false;
  for (it=theVec.begin(); it!=theVec.end(); ++it){
    const jpcRes* jpcIt=(*it).get();
    if (0==jpcIt) ErrMsg(fatal) << "shared object containss 0 pointer!!!" << endmsg;
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


void pbarpStates::print(std::ostream& os) const{
  os << "initital states of the pbar p annihilation for Jmax = " << _jmax << " are: " << std::endl; 
  AbsStates::print(os);

}

