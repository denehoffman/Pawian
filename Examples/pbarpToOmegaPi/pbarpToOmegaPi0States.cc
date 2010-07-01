#include "Examples/pbarpToOmegaPi/pbarpToOmegaPi0States.hh"

#include "Utils/MathUtils.hh"
#include "ErrLogger/ErrLineLog.hh"
#include "PwaUtils/pbarpStates.hh"
 

pbarpToOmegaPi0States::pbarpToOmegaPi0States(boost::shared_ptr<pbarpStates> pbarpStates):
  AbsStates(),
  _omegaJPC(1, -1, -1),
  _piJPC(0, -1, 1),
  _pbarpStatesAll(pbarpStates)
{
  if (0==_pbarpStatesAll) ErrMsg(fatal) << "_pbarpStatesAll==0 !!!" << endmsg;
  calcJPCs();
}


pbarpToOmegaPi0States::~pbarpToOmegaPi0States(){
}

bool pbarpToOmegaPi0States::calcJPCs(){
  

  std::vector< boost::shared_ptr<const jpcRes> > thepbarJPCs=_pbarpStatesAll->jpcStates();
  //  if (0==thepbarJPCs) ErrMsg(fatal) << "_pbarpStatesAll->jpcStates()==0 !!!" << endmsg;

  std::vector< boost::shared_ptr<const jpcRes> >::const_iterator it;

  for ( it=thepbarJPCs.begin(); it!=thepbarJPCs.end(); ++it){
    // first: check C-parity
    if ( (*it)->C != (_omegaJPC.C*_piJPC.C)) continue;

    //second: extract all LS combinations
    vector<LS> motherLSs=GetValidLS((*it)->J, (*it)->P, _omegaJPC.J, _omegaJPC.P, _piJPC.J, _piJPC.P);
    int num_ls = (int) motherLSs.size();
    if (num_ls==0) continue;
    
    ErrMsg(debugging) << "valid LS combinations for JP(pbar p -> omega pi0)=" << (*it)->J <<" " << (*it)->P << endmsg;

    boost::shared_ptr<const jpcRes> jpcPtr=(*it);
    _jpcStates.push_back(jpcPtr);

    for(int ls = 0; ls < num_ls; ls++){
	 Spin L= motherLSs[ls].L;
	 Spin S= motherLSs[ls].S;
	 ErrMsg(debugging) << "L=" << L <<" S=" << S << endmsg;

	 boost::shared_ptr<const JPCLS> tmpJPCLS(new JPCLS(jpcPtr, L, S) );
	 _allJPCLS.push_back(tmpJPCLS);

	 for(Spin M = -S; M <= S; M++){
          double Clebschg=Clebsch(L,0,S,M, (*it)->J,M);
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

