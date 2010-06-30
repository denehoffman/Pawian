#include "Examples/pbarpToOmegaPi/pbarpToOmegaPi0States.hh"

#include "Utils/MathUtils.hh"
#include "ErrLogger/ErrLineLog.hh"
#include "PwaUtils/pbarpStates.hh"
 

pbarpToOmegaPi0States::pbarpToOmegaPi0States(boost::shared_ptr<pbarpStates> pbarpStates):
  AbsStates(),
  _omegaJPC(1, -1, -1),
  _piJPC(0, -1, 1),
  _pbarpStates(pbarpStates)
{
  if (0==_pbarpStates) ErrMsg(fatal) << "_pbarpStates==0 !!!" << endmsg;
  calcJPCs();
}


pbarpToOmegaPi0States::~pbarpToOmegaPi0States(){
}

bool pbarpToOmegaPi0States::calcJPCs(){
  

  std::vector< boost::shared_ptr<const jpcRes> > thepbarJPCs=_pbarpStates->jpcStates();
  //  if (0==thepbarJPCs) ErrMsg(fatal) << "_pbarpStates->jpcStates()==0 !!!" << endmsg;

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

}


void pbarpToOmegaPi0States::print(std::ostream& os) const{
  AbsStates::print(os);

}

