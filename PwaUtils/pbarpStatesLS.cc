#include "PwaUtils/pbarpStatesLS.hh"
#include "Utils/MathUtils.hh"
#include "ErrLogger/ErrLogger.hh"

pbarpStatesLS::pbarpStatesLS():
  AbsStates(),
  _lmax(10),
  _pbarJPC(0.5, -1),
  _pJPC(0.5, 1)
{
   calcStates();
}


pbarpStatesLS::pbarpStatesLS(int lmax):
  AbsStates(),
  _lmax(lmax),
  _pbarJPC(0.5, -1),
  _pJPC(0.5, 1)
{
    calcStates();
}

pbarpStatesLS::~pbarpStatesLS(){
}

bool pbarpStatesLS::calcStates(){
  int jmax=_lmax+1;

  for (int j=0; j<=jmax; j++){
      for (int p=-1; p<=1; p+=2){
	vector<LS> motherLSs=GetValidLS(j, p, _pbarJPC.J,_pbarJPC.P, _pJPC.J, _pJPC.P);
	int num_ls = (int) motherLSs.size();

	for(int ls = 0; ls < num_ls; ls++){

	  Spin L= motherLSs[ls].L;
          if (L>_lmax) continue;

	  Spin S= motherLSs[ls].S;
	  int cparity(pow(-1,int(L+S)));

	  boost::shared_ptr<const jpcRes> jpcPtr(new jpcRes(j,p,cparity));

	  DebugMsg << "L=" << L <<" S=" << S << " ==> C=" << cparity << endmsg;
	  boost::shared_ptr<const JPCLS> jpcLSPtr(new JPCLS(jpcPtr, L, S));

	  for(Spin lampbar = -0.5; lampbar <= 0.5; lampbar++){
	    for(Spin lamp = -0.5; lamp <= 0.5; lamp++){
 
	      double preFactorLS=sqrt(2.*L+1)*Clebsch(L,0,S, lampbar-lamp, j, lampbar-lamp);
	      boost::shared_ptr<const JPCLS> tmpJPCLS( new JPCLS(jpcPtr,L, S, preFactorLS) );              
	      boost::shared_ptr<const JPCLSJJ> jpcLSJJPtr(new JPCLSJJ(tmpJPCLS,_pbarJPC.J, lampbar, _pJPC.J, lamp)); 

 	      if (fabs(jpcLSJJPtr->prefactorAll)>1.e-8){
		_allStates.push_back(jpcLSJJPtr);
		fillVec(jpcPtr, _alljpcRes);
		fillVec(jpcLSPtr, _JPCLS_AllStates);

		if ( fabs(lampbar-0.5) < 1e-8 &&  fabs(lamp-0.5) < 1e-8 ){
		  _lamUpUp_JPCLSJJ_States.push_back(jpcLSJJPtr);
 		  fillVec(tmpJPCLS, _lamUpUp_JPCLS_States);

		  if(jpcLSJJPtr->S==0){
		    fillVec(tmpJPCLS, _singlet_JPCLS_States);
		  }
		  else{ //S==1
		    fillVec(tmpJPCLS, _triplet0_JPCLS_States);
		  }
		}

		else  if ( fabs(lampbar-0.5) < 1e-8 &&  fabs(lamp+0.5) < 1e-8 ){
		  _lamUpDown_JPCLSJJ_States.push_back(jpcLSJJPtr);
		  fillVec(tmpJPCLS, _lamUpDown_JPCLS_States);
		  fillVec(tmpJPCLS, _tripletm1_JPCLS_States);
		}

		else if ( fabs(lampbar+0.5) < 1e-8 &&  fabs(lamp-0.5) < 1e-8 ){
		  _lamDownUp_JPCLSJJ_States.push_back(jpcLSJJPtr);
		  fillVec(tmpJPCLS, _lamDownUp_JPCLS_States);
		  fillVec(tmpJPCLS, _tripletp1_JPCLS_States);
		}

		else if ( fabs(lampbar+0.5) < 1e-8 &&  fabs(lamp+0.5) < 1e-8 ){
		  _lamDownDown_JPCLSJJ_States.push_back(jpcLSJJPtr);

		  fillVec(tmpJPCLS, _lamDownDown_JPCLS_States);

		  if(jpcLSJJPtr->S==0){
		    fillVec(tmpJPCLS, _singlet_JPCLS_States);
		  }
		  else{ //S==1
		    fillVec(tmpJPCLS, _triplet0_JPCLS_States);
		  }

		}

	      } 
	    }
	  }
	}
      }
  }



   return true;
}


void pbarpStatesLS::print(std::ostream& os) const{
  os << "initital states of the pbar p annihilation for Lmax = " << _lmax << " are: " << std::endl; 
  //  AbsStates::print(os);
   std::vector< boost::shared_ptr<const JPCLSJJ > >::const_iterator itJPCLSJJ;
   std::vector< boost::shared_ptr<const JPCLS > >::const_iterator itJPCLS;
   std::vector< boost::shared_ptr<const jpcRes > >::const_iterator itjpcRes;

  os << "\n******** all JPC states ************ " << std::endl;
  for ( itjpcRes=_alljpcRes.begin(); itjpcRes!=_alljpcRes.end(); ++itjpcRes){
    (*itjpcRes)->print(os);
    os << std::endl;
  }

  os << "\n\n******** all JPCLS states ************ " << std::endl;
  for ( itJPCLS=_JPCLS_AllStates.begin(); itJPCLS!=_JPCLS_AllStates.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << std::endl;
  }

  os << "\n\n******** JPCLS states (singlet)  ************ " << std::endl;
  for ( itJPCLS=_singlet_JPCLS_States.begin(); itJPCLS!=_singlet_JPCLS_States.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << std::endl;
  }

  os << "\n******** JPCLS states (triplet 0)  ************ " << std::endl;
  for ( itJPCLS=_triplet0_JPCLS_States.begin(); itJPCLS!=_triplet0_JPCLS_States.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << std::endl;
  }

  os << "\n******** JPCLS states (triplet 1)  ************ " << std::endl;
  for ( itJPCLS=_tripletp1_JPCLS_States.begin(); itJPCLS!=_tripletp1_JPCLS_States.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << std::endl;
  }

  os << "\n******** JPCLS states (triplet -1)  ************ " << std::endl;
  for ( itJPCLS=_tripletm1_JPCLS_States.begin(); itJPCLS!=_tripletm1_JPCLS_States.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << std::endl;
  }

  os << "\n\n******** JPCLS states (lambda(pbar)=1/2 and lambda(p)=1/2)  ************ " << std::endl;
  for ( itJPCLS=_lamUpUp_JPCLS_States.begin(); itJPCLS!=_lamUpUp_JPCLS_States.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << std::endl;
  }

  os << "\n******** JPCLS states (lambda(pbar)=1/2 and lambda(p)=-1/2)  ************ " << std::endl;
  for ( itJPCLS=_lamUpDown_JPCLS_States.begin(); itJPCLS!=_lamUpDown_JPCLS_States.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << std::endl;
  }

  os << "\n******** JPCLS states (lambda(pbar)=-1/2 and lambda(p)=1/2)  ************ " << std::endl;
  for ( itJPCLS=_lamDownUp_JPCLS_States.begin(); itJPCLS!=_lamDownUp_JPCLS_States.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << std::endl;
  }

  os << "\n******** JPCLS states (lambda(pbar)=-1/2 and lambda(p)=-1/2)  ************ " << std::endl;
  for ( itJPCLS=_lamDownDown_JPCLS_States.begin(); itJPCLS!=_lamDownDown_JPCLS_States.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << std::endl;
  }


  os << "\n\n******** all states with lambda(pbar)=1/2 and lambda(p)=1/2************ " << std::endl;
  for ( itJPCLSJJ=_lamUpUp_JPCLSJJ_States.begin(); itJPCLSJJ!=_lamUpUp_JPCLSJJ_States.end(); ++itJPCLSJJ){
    (*itJPCLSJJ)->print(os);
    os << std::endl;
  }

  os << "\n******** all states with lambda(pbar)=1/2 and lambda(p)=-1/2  ************ " << std::endl;
  for ( itJPCLSJJ=_lamUpDown_JPCLSJJ_States.begin(); itJPCLSJJ!=_lamUpDown_JPCLSJJ_States.end(); ++itJPCLSJJ){
    (*itJPCLSJJ)->print(os);
    os << std::endl;
  }

  os << "*\n******* all states with lambda(pbar)=-1/2 and lambda(p)=1/2  ************ " << std::endl;
  for ( itJPCLSJJ=_lamDownUp_JPCLSJJ_States.begin(); itJPCLSJJ!=_lamDownUp_JPCLSJJ_States.end(); ++itJPCLSJJ){
    (*itJPCLSJJ)->print(os);
    os << std::endl;
  }

  os << "\n******** all states with lambda(pbar)=-1/2 and lambda(p)=-1/2  ************ " << std::endl;
  for ( itJPCLSJJ=_lamDownDown_JPCLSJJ_States.begin(); itJPCLSJJ!=_lamDownDown_JPCLSJJ_States.end(); ++itJPCLSJJ){
    (*itJPCLSJJ)->print(os);
    os << std::endl;
  }

}

