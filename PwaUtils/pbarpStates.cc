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
	   
	   bool found=false;
	   std::vector< boost::shared_ptr<const jpcRes> >::const_iterator it;
	   for (it=_jpcStates.begin(); it!=_jpcStates.end(); ++it){
	     const jpcRes* jpcCurrent=jpcPtr.get();
	     const jpcRes* jpcIt=(*it).get();
	     if ( (*jpcIt) == (*jpcCurrent) ){
	       found=true;
	       jpcPtr.reset(); 
	       jpcPtr=(*it);
	       continue;
	     }
	   }
	   if (!found)  _jpcStates.push_back(jpcPtr);
	   
	   boost::shared_ptr<const JPCLS> tmpJPCLS(new JPCLS(jpcPtr, L, S) );
	   _allJPCLS.push_back(tmpJPCLS);
	   
           std::vector<Spin>::const_iterator itM;
	   for ( itM=theMs.begin(); itM!=theMs.end(); ++itM){
	     boost::shared_ptr<const JPCSM> tmpJPCSM( new JPCSM(jpcPtr, S, (*itM) ));
	     _allJPCSM .push_back(tmpJPCSM);
	     
	     double Clebschg1=Clebsch(L,0,S,(*itM), j, (*itM));
	     boost::shared_ptr<const JPCLSM> tmpJPCLSM(new JPCLSM(tmpJPCLS, (*itM), Clebschg1) );
	     _allStates.push_back(tmpJPCLSM);
             if (S==0) _singletStates.push_back(tmpJPCLSM);
	     else if (S==1){
	       if ((*itM)==0) _tripletM0States.push_back(tmpJPCLSM);
               else if ((*itM)==1) _tripletMp1States.push_back(tmpJPCLSM);
               else if ((*itM)==-1) _tripletMm1States.push_back(tmpJPCLSM);
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



void pbarpStates::print(std::ostream& os) const{
  os << "initital states of the pbar p annihilation for Jmax = " << _jmax << " are: " << std::endl; 
  AbsStates::print(os);

}

