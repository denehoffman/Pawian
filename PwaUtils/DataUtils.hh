#ifndef _DataUtils_H
#define _DataUtils_H

#include <boost/shared_ptr.hpp>
#include <string>
#include <sstream>
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

  jpcRes(boost::shared_ptr<const jpcRes> theJPC){
    J=theJPC->J;
    P=theJPC->P;
    C=theJPC->C;
  }

  virtual bool operator==(const jpcRes& compare) const{
    bool result=false;
    if ( fabs(J-compare.J)<1e-8 && P==compare.P && C==compare.C) result=true;
    return result;
  }


  virtual bool operator<(const jpcRes& compare) const{
    bool result=false;
    if ( J < compare.J) result=true;
    else if (J == compare.J){
      if ( P < compare.P) result=true;
      else if (P == compare.P){
	if ( C < compare.C) result=true;
      }
    }
    return result; 
  }

  virtual std::string name() const{
    std::stringstream tmpStrStreamJ;
    tmpStrStreamJ << J;
    std::stringstream tmpStrStreamP;
    tmpStrStreamP << P;   
    std::stringstream tmpStrStreamC;
    tmpStrStreamC << C;   
    std::string result="J"+tmpStrStreamJ.str()+"P"+tmpStrStreamP.str()+"C"+tmpStrStreamC.str();
    return result;
  }

  virtual void print(std::ostream& os) const{
    os <<"J=" << J << "\tP=" << P << "\tC=" << C;   
  }
};


struct JPCLS : public jpcRes{
  Spin L;
  Spin S;
  JPCLS(boost::shared_ptr<const jpcRes> theJPC, const Spin& theL, const Spin& theS): jpcRes(theJPC){
    L=theL;
    S=theS;
  }

  JPCLS(boost::shared_ptr<const JPCLS> theJPCLS): jpcRes(theJPCLS->J, theJPCLS->P, theJPCLS->C){
    L=theJPCLS->L;
    S=theJPCLS->S;
  }

  virtual bool operator==(const jpcRes& compare) const{
    return jpcRes::operator==(compare);
  }

  virtual bool operator==(const JPCLS& compare) const{
    bool result=false;
    if ( fabs(J-compare.J)<1e-8 && P==compare.P && C==compare.C && fabs(L-compare.L)<1e-8 && fabs(S-compare.S)<1e-8 ) result=true;
    return result;
  }

  virtual bool operator<(const JPCLS& compare) const{
    bool result=false;

    if ( J < compare.J) result=true;
    else if (J == compare.J){
      if ( P < compare.P) result=true;
      else if (P == compare.P){
	if ( C < compare.C) result=true;
	else if (C == compare.C){
	  if ( L < compare.L) result=true;
	  else if (L == compare.L){
	    if ( S < compare.S) result=true;
	  }
	}
      }
    }
  
    return result; 
}  

  virtual std::string name() const{
    std::string result=jpcRes::name();
    std::stringstream tmpStrStreamL;
    tmpStrStreamL << L;
    std::stringstream tmpStrStreamS;
    tmpStrStreamS << S;   

    result+="L"+tmpStrStreamL.str()+"S"+tmpStrStreamS.str();
    return result;
  }


  virtual void print(std::ostream& os) const{
    jpcRes::print(os);
    os <<"\tL=" << L << "\tS=" << S;   
  }
};


struct JPCSM : public jpcRes{
  Spin S;
  Spin M;
  JPCSM(boost::shared_ptr<const jpcRes> theJPC, const Spin& theS, const Spin& theM): jpcRes(theJPC){
    S=theS;
    M=theM;
  }

  virtual bool operator==(const jpcRes& compare) const{
    return jpcRes::operator==(compare);
  }

  virtual bool operator==(const JPCSM& compare) const{
    bool result=jpcRes::operator==(compare);;
    if ( fabs(S-compare.S)>1e-8 || fabs(M-compare.M)>1e-8 ) result=false;
    return result;
  }

  virtual bool operator<(const JPCSM& compare) const{
    bool result=false;
    if ( J < compare.J) result=true;
    else if (J == compare.J){
      if ( P < compare.P) result=true;
      else if (P == compare.P){
	if ( C < compare.C) result=true;
	else if (C == compare.C){
	  if ( S < compare.S) result=true;
	  else if (S == compare.S){
	    if ( M < compare.M) result=true;
	  }
	}
      }
    }

    return result; 
  }

  virtual std::string name() const{
    std::string result=jpcRes::name();
    std::stringstream tmpStrStreamS;
    tmpStrStreamS << S;
    std::stringstream tmpStrStreamM;
    tmpStrStreamM << M;   
    
    result+="S"+tmpStrStreamS.str()+"M"+tmpStrStreamM.str();
    return result;
  }
  
  void print(std::ostream& os) const{
    jpcRes::print(os);
    os <<"\tS=" << S << "\tM=" << M
       << std::endl;   
  }
};


struct JPCLSM : public JPCLS{
  boost::shared_ptr<const jpcRes> jpc;
  Spin M;
  double ClebschG;

  JPCLSM(boost::shared_ptr<const JPCLS> theJPCLS, const Spin theM, const double theClebschG) : JPCLS(theJPCLS){
    M=theM;
    ClebschG=theClebschG;
  }

  virtual bool operator==(const jpcRes& compare) const{
    return jpcRes::operator==(compare);
  }

  virtual bool operator==(const JPCLSM& compare) const{
    bool result=jpcRes::operator==(compare);;
    if ( fabs(L-compare.L)>1e-8 || fabs(S-compare.S)>1e-8 || fabs(M-compare.M)>1e-8 ) result=false;
    return result;
  }

  virtual bool operator<(const JPCLSM& compare) const{
    bool result=false;
    if ( J < compare.J) result=true;
    else if (J == compare.J){
      if ( P < compare.P) result=true;
      else if (P == compare.P){
	if ( C < compare.C) result=true;
	else if (C == compare.C){
	  if ( L < compare.L) result=true;
	  else if (L == compare.L){
	    if ( S < compare.S) result=true;
	    else if (S == compare.S){
	      if ( M < compare.M) result=true;
	    }
	  }
	}
      }
    } 
    return result; 
  } 


  virtual std::string name() const{
    std::string result=jpcRes::name();
    std::stringstream tmpStrStreamL;
    tmpStrStreamL << L;
    std::stringstream tmpStrStreamS;
    tmpStrStreamS << S;
    std::stringstream tmpStrStreamM;
    tmpStrStreamS << M;   
    
    result+="L"+tmpStrStreamL.str()+"S"+tmpStrStreamS.str()+"M"+tmpStrStreamM.str();
    return result;
  }


  virtual void print(std::ostream& os) const{
    JPCLS::print(os);
    os <<"\tM=" << M 
       <<"\tClebschGordan=" << ClebschG
       << std::endl;   
  }
};

#endif /* _DataUtils_H */
