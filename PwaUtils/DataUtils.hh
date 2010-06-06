#ifndef _DataUtils_H
#define _DataUtils_H

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

struct jpcRes
{
  Spin J;
  int P;
  int C;
  jpcRes(Spin j=0, int p=1, int c=1){
    J=j;
    P=p;
    C=c;
  }
};

/// L,S,M  combination
struct LSM{
  Spin L;
  Spin S;
  Spin M;
  LSM(int l=0, Spin s=0, int m=0){
    L=l;
    S=s;
    M=m;
  }
};

struct PbarP{
  jpcRes jpc;
  LSM lsm;
  double ClebschG;
  PbarP(const jpcRes& theJPC, const LSM& theLSM, const double theClebschG){
    jpc=theJPC;
    lsm=theLSM;
    ClebschG=theClebschG;
  }

  void print(std::ostream& os) const{
    os <<"J=" << jpc.J << "\tP=" << jpc.P << "\tC=" << jpc.C 
       <<"\tL=" << lsm.L <<"\tS=" << lsm.S <<"\tlambda=" << lsm.M
       <<"\tClebschGordan=" << ClebschG
       << std::endl;   
  }
};

#endif /* _DataUtils_H */
