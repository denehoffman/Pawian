#include "Examples/MATpbarpToOmegaPi/pbarpToOmegaPi0States.hh"

#include "Utils/MathUtils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/AbsStates.hh"
 
boost::shared_ptr<pbarpToOmegaPi0States> pbarpToOmegaPi0States::theStates = boost::shared_ptr<pbarpToOmegaPi0States>();

boost::shared_ptr<pbarpToOmegaPi0States> pbarpToOmegaPi0States::getStates()
{
  return theStates;
}

boost::shared_ptr<pbarpToOmegaPi0States> pbarpToOmegaPi0States::getStates(boost::shared_ptr<pbarpStates> pbarpStates)
{
  if( theStates == 0 ){
    boost::shared_ptr<pbarpToOmegaPi0States> newStates(new pbarpToOmegaPi0States(pbarpStates));
    theStates = newStates;
  }
  return theStates;
}

pbarpToOmegaPi0States::pbarpToOmegaPi0States(boost::shared_ptr<pbarpStates> pbarpStates):
  AbsStates(),
  _omegaJPC(1, -1, -1),
  _piJPC(0, -1, 1),
  _pbarpStatesAll(pbarpStates)
{
  if (0==_pbarpStatesAll) { Alert << "_pbarpStatesAll==0 !!!" << endmsg; exit(1); }
  calcJPCs();
}

pbarpToOmegaPi0States::~pbarpToOmegaPi0States(){
}

vector<LS> pbarpToOmegaPi0States::myGetValidLS(const Spin &__j,int __parity,const Spin &__s1,int __p1,
		      const Spin &__s2,int __p2){

  vector<LS> valid_ls;
  LS ls;

  for(Spin S = abs(__s1 - __s2); S <= (__s1 + __s2); S++){
    for(int L = (int)abs(__j - S); L <= (int)(__j + S); L++){
      // check the parity
      if(abs(__p1*__p2*std::pow(-1.,(int)L) - __parity) < 1.e-5) {
	ls.L = L;
	ls.S = S;
	valid_ls.push_back(ls);
      }
    }
  }
  return valid_ls;
}

bool pbarpToOmegaPi0States::calcJPCs(){
  

  std::vector< boost::shared_ptr<const jpcRes> > thepbarJPCs=_pbarpStatesAll->jpcStates();
  //  if (0==thepbarJPCs) { Alert << "_pbarpStatesAll->jpcStates()==0 !!!" << endmsg; exit(1); }

  std::vector< boost::shared_ptr<const jpcRes> >::const_iterator it;

  for ( it=thepbarJPCs.begin(); it!=thepbarJPCs.end(); ++it){
    // first: check C-parity
    if ( (*it)->C != (_omegaJPC.C*_piJPC.C)) continue;

    //second: extract all LS combinations
    vector<LS> motherLSs=myGetValidLS((*it)->J, (*it)->P, _omegaJPC.J, _omegaJPC.P, _piJPC.J, _piJPC.P);
    int num_ls = (int) motherLSs.size();
    if (num_ls==0) continue;
    
    DebugMsg << "valid LS combinations for JP(pbar p -> omega pi0)=" << (*it)->J <<" " << (*it)->P << endmsg;

    boost::shared_ptr<const jpcRes> jpcPtr=(*it);
    _jpcStates.push_back(jpcPtr);

    for(int ls = 0; ls < num_ls; ls++){
	 Spin L= motherLSs[ls].L;
	 Spin S= motherLSs[ls].S;
	 DebugMsg << "L=" << L <<" S=" << S << endmsg;

	 boost::shared_ptr<const JPCLS> tmpJPCLS(new JPCLS(jpcPtr, L, S) );
	 _allJPCLS.push_back(tmpJPCLS);

	 for(Spin M = -S; M <= S; M++){
          Spin myJ((*it)->J.Numerator(), (*it)->J.Denominator());
          //double Clebschg=Clebsch(L,0,S,M, (*it)->J,M);
          double Clebschg=Clebsch(L,0,S,M, myJ,M);
	  boost::shared_ptr<const JPCLSM> tmpJPCLSM(new JPCLSM(tmpJPCLS, M, Clebschg) );
          _allStates.push_back(tmpJPCLSM);   
	 }
    }
  }

  _pbarpSingletToOmegaPi=_pbarpStatesAll->extractSingletStates(_jpcStates);
  _pbarpTripletM0ToOmegaPi=_pbarpStatesAll->extractTripletM0States(_jpcStates);
  _pbarpTripletM1ToOmegaPi=_pbarpStatesAll->extractTripletM1States(_jpcStates);
  _JPCLSomegaProdSinglet=extractJPCLSStates(_pbarpSingletToOmegaPi, _allJPCLS);
  _JPCLSomegaProdTripletM0=extractJPCLSStates(_pbarpTripletM0ToOmegaPi, _allJPCLS);
  _JPCLSomegaProdTripletM1=extractJPCLSStates(_pbarpTripletM1ToOmegaPi, _allJPCLS);
}

std::vector< boost::shared_ptr<const JPCLS> > pbarpToOmegaPi0States::extractJPCLSStates(std::vector< boost::shared_ptr<const jpcRes> > theJPCRequests, std::vector< boost::shared_ptr<const JPCLS> > jpclsStatesAll) const{

  std::vector< boost::shared_ptr<const JPCLS> > result;

  std::vector< boost::shared_ptr<const jpcRes> >::const_iterator itJPC;
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itAllStates;
  for ( itJPC=theJPCRequests.begin(); itJPC!=theJPCRequests.end(); ++itJPC){
    const jpcRes* jpcRequest=(*itJPC).get();

    for ( itAllStates=jpclsStatesAll.begin(); itAllStates!=jpclsStatesAll.end(); ++itAllStates){
      const JPCLS* jpclsCurrent=(*itAllStates).get();
      if (( *jpclsCurrent) ==  (*jpcRequest)){
        result.push_back(*itAllStates);
        continue;
      } 
    }
  }
  return result;
}

void pbarpToOmegaPi0States::print(std::ostream& os) const{
  os << "pbarpToOmegaPi0States::print\n"
     << "The contributed JPC initial states are:" << std::endl;
  
  std::vector< boost::shared_ptr<const jpcRes > >::const_iterator itJPC;
  for ( itJPC=_jpcStates.begin(); itJPC!=_jpcStates.end(); ++itJPC){
    (*itJPC)->print(os);
    os << std::endl;
  }

  os << "The initial JPC singlet states are:" << std::endl;
  for ( itJPC=_pbarpSingletToOmegaPi.begin(); itJPC!=_pbarpSingletToOmegaPi.end(); ++itJPC){
    (*itJPC)->print(os);
    os << std::endl;
  }

  os << "The initial JPC triplet M=0 states are:" << std::endl;
  for ( itJPC=_pbarpTripletM0ToOmegaPi.begin(); itJPC!=_pbarpTripletM0ToOmegaPi.end(); ++itJPC){
    (*itJPC)->print(os);
    os << std::endl;
  }

  os << "The initial JPC triplet M=1 states are:" << std::endl;
  for ( itJPC=_pbarpTripletM1ToOmegaPi.begin(); itJPC!=_pbarpTripletM1ToOmegaPi.end(); ++itJPC){
    (*itJPC)->print(os);
    os << std::endl;
  }

  os << "The contributed JPC LS combinations for the omega production are:" << std::endl;
  std::vector< boost::shared_ptr<const JPCLS > >::const_iterator itJPCLS;
  for ( itJPCLS=_allJPCLS.begin(); itJPCLS!=_allJPCLS.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << std::endl;
  }

  os << "The fit parameter for the pbarp singet states are:" << std::endl;
  for ( itJPCLS=_JPCLSomegaProdSinglet.begin(); itJPCLS!=_JPCLSomegaProdSinglet.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << std::endl;
  }

  os << "The fit parameter for the pbarp triplet M=0 states are:" << std::endl;
  for ( itJPCLS=_JPCLSomegaProdTripletM0.begin(); itJPCLS!=_JPCLSomegaProdTripletM0.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << std::endl;
  }

  os << "The fit parameter for the pbarp triplet M=+-1 states are:" << std::endl;
  for ( itJPCLS=_JPCLSomegaProdTripletM1.begin(); itJPCLS!=_JPCLSomegaProdTripletM1.end(); ++itJPCLS){
    (*itJPCLS)->print(os);
    os << std::endl;
  }
}



