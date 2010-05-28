#include "PwaUtils/pbarpStates.hh"
#include "ErrLogger/ErrLineLog.hh"


 

pbarpStates::pbarpStates():
  _jmax(10),
  _pbarJPC(0.5, -1),
  _pJPC(0.5, 1){
  calcJPCs();
}


pbarpStates::pbarpStates(int jmax):
  _jmax(jmax),
  _pbarJPC(0.5, -1),
  _pJPC(0.5, 1)
{
  calcJPCs();
}

pbarpStates::~pbarpStates(){
  std::vector<PbarP*>::iterator it;
  for ( it=_theStates.begin(); it!=_theStates.end(); ++it){
    if (0!= (*it)){
      delete (*it);
      (*it)=0;
    }
  }
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
	 int cparity(pow(-1,L+S));
         ErrMsg(debugging) << "L=" << L <<" S=" << S << " ==> C=" << cparity << endmsg;
         for(Spin M = -S; M <= S; M++){
           double Clebschg=Clebsch(L,0,S,M, j,M);
	   ErrMsg(debugging) << "Clebsch(L,0,S,M=" << M << ", j,M=" << M << "): " << Clebschg << endmsg;
           if (fabs(Clebschg)>1e-8){
             jpcRes theJPC(j,p,cparity);
             LSM theLSM(L,S,M);
	     _theStates.push_back(new PbarP(theJPC, theLSM, Clebschg));
	   }
	 }
       }
       ErrMsg(debugging) << "\n" << endmsg;
      }
    }
}

void pbarpStates::print(std::ostream& os) const{
  os << "initital states of the pbar p annihilation for Jmax = " << _jmax << " are: " << std::endl; 

  std::vector<PbarP*>::const_iterator it;
   for ( it=_theStates.begin(); it!=_theStates.end(); ++it){
    if (0!= (*it)){
	os <<"J=" << (*it)->jpc.J << " P=" << (*it)->jpc.P << " C=" << (*it)->jpc.C 
	   <<" L=" << (*it)->lsm.L <<" S=" << (*it)->lsm.S <<" lambda=" << (*it)->lsm.M
	   <<" ClebschGordan=" << (*it)->ClebschG 
	   << std::endl;
      
    }
   }
}
