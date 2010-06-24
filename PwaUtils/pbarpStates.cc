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
         for(Spin M = -S; M <= S; M++){
           double Clebschg=Clebsch(L,0,S,M, j,M);
	   ErrMsg(debugging) << "Clebsch(L,0,S,M=" << M << ", j,M=" << M << "): " << Clebschg << endmsg;
           if (fabs(Clebschg)>1e-8){

 	     boost::shared_ptr<const jpcRes> jpcPtr(new jpcRes(j,p,cparity));

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

             boost::shared_ptr<const JPCSM> tmpJPCSM(new JPCSM(jpcPtr, S, M) );
	     found=false;
	     std::vector< boost::shared_ptr<const JPCSM> >::const_iterator it1;
	     for (it1=_allJPCSM.begin(); it1!=_allJPCSM.end(); ++it1){
	       const JPCSM* jpcsmCurrent=tmpJPCSM.get();
               const JPCSM* jpcsmIt=(*it1).get();
	       if ( (*jpcsmIt) == (*jpcsmCurrent) ){
		 found=true;
		 tmpJPCSM.reset(); 
                 tmpJPCSM=(*it1);
		 continue;
	       }
	     }
             if (!found)  _allJPCSM .push_back(tmpJPCSM);

             LSM theLSM(L,S,M);

             boost::shared_ptr<const JPCSML> tmpJPCSML(new JPCSML(jpcPtr, theLSM, Clebschg) );

	     //and now fill the vectors
	     _allStates.push_back(tmpJPCSML);
             if (S==0) _singletStates.push_back(tmpJPCSML);
	     else if (S==1){
	       if (M==0) _tripletM0States.push_back(tmpJPCSML);
               else if (M==1) _tripletMp1States.push_back(tmpJPCSML);
               else if (M==-1) _tripletMm1States.push_back(tmpJPCSML);
	       else ErrMsg(fatal) << "pbar p state with S=" << S << " and M="<< M <<" cannot exitst!!!" << cparity << endmsg;
	     }
	     else ErrMsg(fatal) << "pbar p state with S=" << S << " cannot exitst!!!" << cparity << endmsg;
	   }
	 }
       }
       ErrMsg(debugging) << "\n" << endmsg;
      }
    }
}


void pbarpStates::print(std::ostream& os) const{
  os << "initital states of the pbar p annihilation for Jmax = " << _jmax << " are: " << std::endl; 
  AbsStates::print(os);

}

