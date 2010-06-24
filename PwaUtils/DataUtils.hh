#ifndef _DataUtils_H
#define _DataUtils_H

#include <boost/shared_ptr.hpp>
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

  bool operator==(const jpcRes& compare) const{
    bool result=false;
    if ( fabs(J-compare.J)<1e-8 && P==compare.P && C==compare.C) result=true;
    return result;
  }

  bool operator<(const jpcRes& compare){
    if ( J < compare.J) return true;
    else if ( P < compare.P) return true;
    else if ( C < compare.C) return true;
    return false; 
  }

  void print(std::ostream& os) const{
    os <<"J=" << J << "\tP=" << P << "\tC=" << C;   
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

  void print(std::ostream& os) const{
    os <<"L=" << L << "\tS=" << S << "\tM=" << M;   
  }
};


struct JPCSM{
  boost::shared_ptr<const jpcRes> jpc;
  Spin S;
  Spin M;
  JPCSM(boost::shared_ptr<const jpcRes> theJPC, const Spin& theS, const Spin& theM){
    jpc=theJPC;
    S=theS;
    M=theM;
  }

  bool operator==(const JPCSM& compare) const{
    bool result=false;
    const jpcRes* jpcRequest=compare.jpc.get();
    const jpcRes* jpcThis=jpc.get();
    if ( (*jpcThis) == (*jpcRequest) && fabs(S-compare.S)<1e-8 && fabs(M-compare.M)<1e-8 ) result=true;
    return result;
  }

  void print(std::ostream& os) const{
    jpc->print(os);
    os <<"\tS=" << S << "\tM=" << M
       << std::endl;   
  }
};


struct JPCSML{
  boost::shared_ptr<const jpcRes> jpc;
  LSM lsm;
  double ClebschG;
  JPCSML(boost::shared_ptr<const jpcRes> theJPC, const LSM& theLSM, const double theClebschG){
    jpc=theJPC;
    lsm=theLSM;
    ClebschG=theClebschG;
  }

  void print(std::ostream& os) const{
    jpc->print(os);
    os <<"\t";
    lsm.print(os);
    os <<"\tClebschGordan=" << ClebschG
       << std::endl;   
  }
};

#endif /* _DataUtils_H */
