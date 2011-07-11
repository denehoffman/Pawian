#include "Examples/pbarpToOmegaPiLS/pbarpToOmegaPi0StatesLS.hh"

#include "Utils/MathUtils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/pbarpStates.hh"
 

pbarpToOmegaPi0StatesLS::pbarpToOmegaPi0StatesLS(int Lmax):
   pbarpStatesLS(Lmax),
  _omegaJPC(1, -1, -1),
  _piJPC(0, -1, 1)
{
//   _lmax=Lmax;
  _alljpcRes.clear();
  _allStates.clear();
  _lamUpUp_JPCLSJJ_States.clear();
  _lamUpDown_JPCLSJJ_States.clear();
  _lamDownUp_JPCLSJJ_States.clear();
  _lamDownDown_JPCLSJJ_States.clear();
  _lamUpUp_JPCLS_States.clear();
  _lamUpDown_JPCLS_States.clear();
  _lamDownUp_JPCLS_States.clear();
  _lamDownDown_JPCLS_States.clear();
  _JPCLS_AllStates.clear();
  _singlet_JPCLS_States.clear();
  _triplet0_JPCLS_States.clear();
  _tripletp1_JPCLS_States.clear();
  _tripletm1_JPCLS_States.clear();
   calcStates();
}



pbarpToOmegaPi0StatesLS::~pbarpToOmegaPi0StatesLS(){
}

bool pbarpToOmegaPi0StatesLS::calcStates(){
  
  pbarpStatesLS allpbarpStateLS(_lmax);

  std::vector< boost::shared_ptr<const jpcRes> > thepbarpJPCs=allpbarpStateLS.jpcStates();

//   if ( 0 == thepbarpJPCs ) { Alert << "thepbarpJPCs->jpcStates()==0 !!!" << endmsg; exit(1); }  

  std::vector< boost::shared_ptr<const jpcRes> >::const_iterator it;

  for ( it=thepbarpJPCs.begin(); it!=thepbarpJPCs.end(); ++it){
    // first: check C-parity
    if ( (*it)->C != (_omegaJPC.C*_piJPC.C)) continue;

    //second: extract all LS combinations
    vector<LS> motherLSs=GetValidLS((*it)->J, (*it)->P, _omegaJPC.J, _omegaJPC.P, _piJPC.J, _piJPC.P);
    int num_ls = (int) motherLSs.size();
    if (num_ls==0) continue;

    DebugMsg << "valid LS combinations for JP(pbar p -> omega pi0)=" << (*it)->J <<" " << (*it)->P << endmsg;

    boost::shared_ptr<const jpcRes> jpcPtr=(*it);

    for(int ls = 0; ls < num_ls; ls++){
      Spin l= motherLSs[ls].L;
      Spin s= motherLSs[ls].S;
      DebugMsg << "l=" << l <<" s=" << s << endmsg;
      
      boost::shared_ptr<const JPCLS> tmpJPCLS(new JPCLS(jpcPtr, l, s) );
      
      for(Spin M = -s; M <= s; M++){
	double Clebschg=Clebsch(l,0,s,M, (*it)->J, M);
        if ( fabs(Clebschg)>1e-8){
	  fillVec(jpcPtr, _alljpcRes);
	  fillVec(tmpJPCLS, _allOmegaPiStates);

	}
      }
    }
  }

  // now fill all pbarp states vectors
  std::vector< boost::shared_ptr<const JPCLS> > allJPCLSStates=allpbarpStateLS.all_JPCLS_States();
  _JPCLS_AllStates= extractStates(_alljpcRes, allJPCLSStates);

  std::vector< boost::shared_ptr<const JPCLS> > uuJPCLSStates=allpbarpStateLS.lamUpUp_JPCLS_States();
  _lamUpUp_JPCLS_States= extractStates(_alljpcRes, uuJPCLSStates);

  std::vector< boost::shared_ptr<const JPCLS> > udJPCLSStates=allpbarpStateLS.lamUpDown_JPCLS_States();
  _lamUpDown_JPCLS_States= extractStates(_alljpcRes, udJPCLSStates);

  std::vector< boost::shared_ptr<const JPCLS> > duJPCLSStates=allpbarpStateLS.lamDownUp_JPCLS_States();
  _lamDownUp_JPCLS_States= extractStates(_alljpcRes, duJPCLSStates);

   std::vector< boost::shared_ptr<const JPCLS> > ddJPCLSStates=allpbarpStateLS.lamDownDown_JPCLS_States();
  _lamDownDown_JPCLS_States= extractStates(_alljpcRes, ddJPCLSStates);

  std::vector< boost::shared_ptr<const JPCLSJJ> > allJPCLSJJStates=allpbarpStateLS.all_JPCLSJJ_States();  
  _allStates= extractStates(_alljpcRes, allJPCLSJJStates);

  std::vector< boost::shared_ptr<const JPCLSJJ> > uuJPCLSJJStates=allpbarpStateLS.lamUpUp_JPCLSJJ_States();
  _lamUpUp_JPCLSJJ_States= extractStates(_alljpcRes, uuJPCLSJJStates);

  std::vector< boost::shared_ptr<const JPCLSJJ> > udJPCLSJJStates=allpbarpStateLS.lamUpDown_JPCLSJJ_States();
  _lamUpDown_JPCLSJJ_States= extractStates(_alljpcRes, udJPCLSJJStates);

  std::vector< boost::shared_ptr<const JPCLSJJ> > duJPCLSJJStates=allpbarpStateLS.lamDownUp_JPCLSJJ_States();
  _lamDownUp_JPCLSJJ_States= extractStates(_alljpcRes, duJPCLSJJStates);

  std::vector< boost::shared_ptr<const JPCLSJJ> > ddJPCLSJJStates=allpbarpStateLS.lamDownDown_JPCLSJJ_States();
  _lamDownDown_JPCLSJJ_States= extractStates(_alljpcRes, ddJPCLSJJStates);  


  std::vector< boost::shared_ptr<const JPCLS> > singletJPCLS=allpbarpStateLS.singlet_JPCLS_States();
  _singlet_JPCLS_States= extractStates(_alljpcRes, singletJPCLS);

  std::vector< boost::shared_ptr<const JPCLS> > triplet0JPCLS=allpbarpStateLS.triplet0_JPCLS_States();
  _triplet0_JPCLS_States= extractStates(_alljpcRes, triplet0JPCLS);

  std::vector< boost::shared_ptr<const JPCLS> > tripletp1JPCLS=allpbarpStateLS.tripletp1_JPCLS_States();
  _tripletp1_JPCLS_States= extractStates(_alljpcRes, tripletp1JPCLS);

  std::vector< boost::shared_ptr<const JPCLS> > tripletm1JPCLS=allpbarpStateLS.tripletm1_JPCLS_States();
  _tripletm1_JPCLS_States= extractStates(_alljpcRes, tripletm1JPCLS);

  // now fill all initial+decay states 
  
  fillStatesInitialDec( _JPCLS_AllStates, _allOmegaPiStates, _allJPCLSlsStates);  
  fillStatesInitialDec( _singlet_JPCLS_States, _allOmegaPiStates,_singlet_JPCLS_StatesDec);
  fillStatesInitialDec( _triplet0_JPCLS_States, _allOmegaPiStates,_triplet0_JPCLS_StatesDec);
  fillStatesInitialDec( _tripletp1_JPCLS_States, _allOmegaPiStates,_tripletp1_JPCLS_StatesDec);
  fillStatesInitialDec( _tripletm1_JPCLS_States, _allOmegaPiStates,_tripletm1_JPCLS_StatesDec);

  fillStatesInitialDec( _lamUpUp_JPCLSJJ_States, _allOmegaPiStates, _lamUpUp_JPCLSJJ_StatesDec);
  fillStatesInitialDec( _lamUpDown_JPCLSJJ_States, _allOmegaPiStates, _lamUpDown_JPCLSJJ_StatesDec);
  fillStatesInitialDec( _lamDownUp_JPCLSJJ_States, _allOmegaPiStates, _lamDownUp_JPCLSJJ_StatesDec);
  fillStatesInitialDec( _lamDownDown_JPCLSJJ_States, _allOmegaPiStates, _lamDownDown_JPCLSJJ_StatesDec);
  return true;
}

void pbarpToOmegaPi0StatesLS::print(std::ostream& os) const{

   pbarpStatesLS::print(os);

  os << "\n \n \n ************* pbar p -> omega pi0 states for Lmax = " << _lmax << " are: *************" << std::endl;

   std::vector< boost::shared_ptr<const JPCLS > >::const_iterator itJPCLS;

  os << "\n\n******* all JPCls states ************ " << std::endl;
  for ( itJPCLS=_allOmegaPiStates.begin(); itJPCLS!=_allOmegaPiStates.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << std::endl;
  }

   std::vector< boost::shared_ptr<const JPCLSls > >::const_iterator itJPCLSls;
  os << "\n\n******** all JPCLSls states ************ " << std::endl;
  for ( itJPCLSls=_allJPCLSlsStates.begin(); itJPCLSls!=_allJPCLSlsStates.end(); ++itJPCLSls){
    (*itJPCLSls)->print(os);
    os << std::endl;
  } 

  os << "\n\n******** JPCLSls states (singlet)************ " << std::endl;
  for ( itJPCLSls=_singlet_JPCLS_StatesDec.begin(); itJPCLSls!=_singlet_JPCLS_StatesDec.end(); ++itJPCLSls){
    (*itJPCLSls)->print(os);
    os << std::endl;
  } 

  os << "\n******** JPCLSls states (triplet 0)************ " << std::endl;
  for ( itJPCLSls=_triplet0_JPCLS_StatesDec.begin(); itJPCLSls!=_triplet0_JPCLS_StatesDec.end(); ++itJPCLSls){
    (*itJPCLSls)->print(os);
    os << std::endl;
  } 

  os << "\n******** JPCLSls states (triplet 1)************ " << std::endl;
  for ( itJPCLSls=_tripletp1_JPCLS_StatesDec.begin(); itJPCLSls!=_tripletp1_JPCLS_StatesDec.end(); ++itJPCLSls){
    (*itJPCLSls)->print(os);
    os << std::endl;
  }

  os << "\n******** JPCLSls states (triplet -1)************ " << std::endl;
  for ( itJPCLSls=_tripletm1_JPCLS_StatesDec.begin(); itJPCLSls!=_tripletm1_JPCLS_StatesDec.end(); ++itJPCLSls){
    (*itJPCLSls)->print(os);
    os << std::endl;
  }

  std::vector< boost::shared_ptr<const JPCLSJJls > >::const_iterator itJPCLSJJls;

  os << "\n******** JPCLSJJls states (lambda(pbar)=1/2 and lambda(p)=1/2)  ************ " << std::endl;
  for ( itJPCLSJJls=_lamUpUp_JPCLSJJ_StatesDec.begin(); itJPCLSJJls!=_lamUpUp_JPCLSJJ_StatesDec.end(); ++itJPCLSJJls){
    (*itJPCLSJJls)->print(os);
    os << std::endl;
  }

  os << "\n******** LS states (lambda(pbar)=1/2 and lambda(p)=-1/2)  ************ " << std::endl;
  for ( itJPCLSJJls=_lamUpDown_JPCLSJJ_StatesDec.begin(); itJPCLSJJls!=_lamUpDown_JPCLSJJ_StatesDec.end(); ++itJPCLSJJls){
    (*itJPCLSJJls)->print(os);
    os << std::endl;
  }

  os << "\n******** LS states (lambda(pbar)=-1/2 and lambda(p)=1/2)  ************ " << std::endl;
  for ( itJPCLSJJls=_lamDownUp_JPCLSJJ_StatesDec.begin(); itJPCLSJJls!=_lamDownUp_JPCLSJJ_StatesDec.end(); ++itJPCLSJJls){
    (*itJPCLSJJls)->print(os);
    os << std::endl;
  }

  os << "\n******** LS states (lambda(pbar)=-1/2 and lambda(p)=1/2)  ************ " << std::endl;
  for ( itJPCLSJJls=_lamDownDown_JPCLSJJ_StatesDec.begin(); itJPCLSJJls!=_lamDownDown_JPCLSJJ_StatesDec.end(); ++itJPCLSJJls){
    (*itJPCLSJJls)->print(os);
    os << std::endl;
  }
}

