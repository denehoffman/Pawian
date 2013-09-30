//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
//									  //
//  This file is part of Pawian.					  //
//									  //
//  Pawian is free software: you can redistribute it and/or modify	  //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	  //
//									  //
//  Pawian is distributed in the hope that it will be useful,		  //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  //
//  GNU General Public License for more details.	      		  //
//									  //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	  //
//									  //
//************************************************************************//

// DataUtils class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"


struct LHData{
   double weightSum;
   long num_mc;
   double logLH_data;
   double LH_mc;

  LHData(){
    weightSum=logLH_data=LH_mc=num_mc=0.0;
  }
};



struct jpcRes {
  Spin J;
  int P;
  int C;

  jpcRes(Spin j=0, int p=1, int c=1) {
    J=j;
    P=p;
    C=c;
  }

  jpcRes(std::shared_ptr<const jpcRes> theJPC) {
    J=theJPC->J;
    P=theJPC->P;
    C=theJPC->C;
  }

  virtual bool operator==(const jpcRes& compare) const {
    bool result=false;
    if ( fabs(J-compare.J)<1e-8 && P==compare.P && C==compare.C) result=true;
    return result;
  }

  virtual bool operator<(const jpcRes& compare) const {
    bool result=false;
    if ( J < compare.J) result=true;
    else if (J == compare.J) {
      if ( P < compare.P) result=true;
      else if (P == compare.P) {
	if ( C < compare.C) result=true;
      }
    }
    return result;
  }
  virtual std::string jpcname() const {
    std::stringstream tmpStrStreamJ;
    tmpStrStreamJ << J;
    std::stringstream tmpStrStreamP;
    tmpStrStreamP << P;

    std::stringstream tmpStrStreamC;
    tmpStrStreamC << C;
    std::string result="J"+tmpStrStreamJ.str()+"P"+tmpStrStreamP.str()+"C"+tmpStrStreamC.str();
    return result;
  }

  virtual std::string name() const {
    return jpcname();
    // std::stringstream tmpStrStreamJ;
    // tmpStrStreamJ << J;
    // std::stringstream tmpStrStreamP;
    // tmpStrStreamP << P;

    // std::stringstream tmpStrStreamC;
    // tmpStrStreamC << C;
    // std::string result="J"+tmpStrStreamJ.str()+"P"+tmpStrStreamP.str()+"C"+tmpStrStreamC.str();
    // return result;
  }

  virtual void print(std::ostream& os) const {
    os <<"JPC=" << J;
    if (P<0.) os << "-";
    else if (P>0.) os <<  "+";
    else os << "0";

    if (C<0.) os << "-";
    else if (C>0.) os <<  "+";
    else os << "0";
  }

  virtual ~jpcRes(){};
};




struct IGJPC : public jpcRes {
  Spin I;
  int G;
  double parityFactor;

  IGJPC(const Spin j, const int p, const int c, const Spin& theI,
	    int theG, const double theParityFactor=0.): jpcRes(j, p, c) {
    I=theI;
    G=theG;
    parityFactor=theParityFactor;
   }

  IGJPC(std::shared_ptr<const jpcRes> theJPC, const Spin& theI,
	    int theG, const double theParityFactor=0.): jpcRes(theJPC) {
    I=theI;
    G=theG;
    parityFactor=theParityFactor;
  }



  IGJPC(std::shared_ptr<const IGJPC> theIGJPC):
    jpcRes(theIGJPC->J, theIGJPC->P, theIGJPC->C) {
    I=theIGJPC->I;
    G=theIGJPC->G;
    parityFactor=theIGJPC->parityFactor;
  }

  virtual bool operator==(const jpcRes& compare) const {
    return jpcRes::operator==(compare);
  }

  virtual bool operator==(const IGJPC& compare) const {
    bool result=false;
    if ( fabs(J-compare.J)<1e-8 && P==compare.P && C==compare.C
	 && fabs(I-compare.I)<1e-8 && fabs(G-compare.G)<1e-8)
      result=true;
    return result;
  }

  virtual bool operator<(const IGJPC& compare) const {
    bool result=false;

    if ( J < compare.J) result=true;
    else if (J == compare.J){
      if ( P < compare.P) result=true;
      else if (P == compare.P){
	if ( C < compare.C) result=true;
	else if (C == compare.C){
	  if ( I < compare.I) result=true;
	  else if (I == compare.I){
	    if ( G < compare.G) result=true;
	  }
	}
      }
    }

    return result;
}

  virtual std::string name() const {
    //    std::string result=jpcRes::name();
    std::stringstream tmpStrStreamI;
    tmpStrStreamI << I;
    std::stringstream tmpStrStreamG;
    tmpStrStreamG << G;

    std::string result="I"+tmpStrStreamI.str()+"G"+tmpStrStreamG.str()+jpcRes::name();
    return result;
  }


  virtual void print(std::ostream& os) const {
    os <<"\tI,G=" << I << "," << G;
    jpcRes::print(os);
    os << "\tparityFactor=" << parityFactor;
  }
};


struct JPClamlam : public jpcRes {
  Spin lam1;
  Spin lam2;
  double parityFactor;

  JPClamlam(const Spin j, const int p, const int c, const Spin& theLam1,
	    const Spin& theLam2, const double theParityFactor=0.): jpcRes(j, p, c) {
    lam1=theLam1;
    lam2=theLam2;
    parityFactor=theParityFactor;
   }

  JPClamlam(std::shared_ptr<const jpcRes> theJPC, const Spin& theLam1,
	    const Spin& theLam2, const double theParityFactor=0.): jpcRes(theJPC) {
    lam1=theLam1;
    lam2=theLam2;
    parityFactor=theParityFactor;
  }

  JPClamlam(std::shared_ptr<const JPClamlam> theJPClamlam):
    jpcRes(theJPClamlam->J, theJPClamlam->P, theJPClamlam->C) {
    lam1=theJPClamlam->lam1;
    lam2=theJPClamlam->lam2;
    parityFactor=theJPClamlam->parityFactor;
  }

  virtual bool operator==(const jpcRes& compare) const {
    return jpcRes::operator==(compare);
  }

  virtual bool operator==(const JPClamlam& compare) const {
    bool result=false;
    if ( fabs(J-compare.J)<1e-8 && P==compare.P && C==compare.C
	 && fabs(lam1-compare.lam1)<1e-8 && fabs(lam2-compare.lam2)<1e-8)
      result=true;
    return result;
  }

  virtual bool operator<(const JPClamlam& compare) const {
    bool result=false;

    if ( J < compare.J) result=true;
    else if (J == compare.J){
      if ( P < compare.P) result=true;
      else if (P == compare.P){
	if ( C < compare.C) result=true;
	else if (C == compare.C){
	  if ( lam1 < compare.lam1) result=true;
	  else if (lam1 == compare.lam1){
	    if ( lam2 < compare.lam2) result=true;
	  }
	}
      }
    }

    return result;
}

  virtual std::string name() const {
    std::string result=jpcRes::name();
    std::stringstream tmpStrStreamLam1;
    tmpStrStreamLam1 << lam1;
    std::stringstream tmpStrStreamLam2;
    tmpStrStreamLam2 << lam2;

    result+="Lama"+tmpStrStreamLam1.str()+"Lamb"+tmpStrStreamLam2.str();
    return result;
  }


  virtual void print(std::ostream& os) const {
    jpcRes::print(os);
    os <<"\tLam1,Lam2=" << lam1 << "," << lam2 <<"\tparityFactor=" << parityFactor;
  }
};


struct JPCLS : public jpcRes {
  Spin L;
  Spin S;
  double preFactor;


  JPCLS(const Spin j, const int p, const int c, const Spin& theL, const Spin& theS,
	const double thePreFactor=0.): jpcRes(j, p, c) {
    L=theL;
    S=theS;
    preFactor=thePreFactor;
  }

  JPCLS(std::shared_ptr<const jpcRes> theJPC, const Spin& theL, const Spin& theS,
	const double thePreFactor=0.): jpcRes(theJPC) {
    L=theL;
    S=theS;
    preFactor=thePreFactor;
  }

  JPCLS(std::shared_ptr<const JPCLS> theJPCLS):
    jpcRes(theJPCLS->J, theJPCLS->P, theJPCLS->C) {
    L=theJPCLS->L;
    S=theJPCLS->S;
    preFactor=theJPCLS->preFactor;
  }

  virtual bool operator==(const jpcRes& compare) const {
    return jpcRes::operator==(compare);
  }

  virtual bool operator==(const JPCLS& compare) const {
    bool result=false;
    if ( fabs(J-compare.J)<1e-8 && P==compare.P && C==compare.C
	 && fabs(L-compare.L)<1e-8 && fabs(S-compare.S)<1e-8)
      result=true;
    return result;
  }

  virtual bool operator<(const JPCLS& compare) const {
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

  virtual std::string name() const {
    std::string result=jpcRes::name();
    std::stringstream tmpStrStreamL;
    tmpStrStreamL << L;
    std::stringstream tmpStrStreamS;
    tmpStrStreamS << S;

    result+="L"+tmpStrStreamL.str()+"S"+tmpStrStreamS.str();
    return result;
  }


  virtual void print(std::ostream& os) const {
    jpcRes::print(os);
    os <<"\tL,S=" << L << "," << S << "\tpreFactorLS=" << preFactor;
  }
};


struct JPCLSJJ : public JPCLS {
  Spin J1;
  Spin Lambda1;

  Spin J2;
  Spin Lambda2;

  double CGLS;
  double CGJJ;
  double prefactorAll;

  JPCLSJJ(std::shared_ptr<const JPCLS> theJPCLS, const Spin& theJ1,
	  const Spin& theLambda1, const Spin& theJ2, const Spin& theLambda2):
    JPCLS(theJPCLS) {
    J1=theJ1;
    Lambda1=theLambda1;
    J2=theJ2;
    Lambda2=theLambda2;
    CGLS=Clebsch(L,0,S, Lambda1-Lambda2, J, Lambda1-Lambda2);
    CGJJ=Clebsch(J1, Lambda1, J2, -Lambda2, S, Lambda1-Lambda2);
    prefactorAll = sqrt( (2.*L+1)/(2.*J+1)  ) * CGLS * CGJJ;
  }

  JPCLSJJ(std::shared_ptr<const JPCLSJJ> theJPCLSJJ): JPCLS(theJPCLSJJ) {
    J1=theJPCLSJJ->J1;
    Lambda1=theJPCLSJJ->Lambda1;
    J2=theJPCLSJJ->J2;
    Lambda2=theJPCLSJJ->Lambda2;
    CGLS=theJPCLSJJ->CGLS;
    CGJJ=theJPCLSJJ->CGJJ;
    prefactorAll=theJPCLSJJ->prefactorAll;
  }

  virtual bool operator==(const jpcRes& compare) const {
    return jpcRes::operator==(compare);
  }

  virtual bool operator==(const JPCLSJJ& compare) const {
    bool result = JPCLS::operator==(compare);
    if ( fabs(J1-compare.J1)>1e-8 || fabs(Lambda1-compare.Lambda1)>1e-8
	 || fabs(J2-compare.J2)>1e-8 || fabs(Lambda2-compare.Lambda2)>1e-8 )
      result=false;
    return result;
  }

  virtual bool operator<(const JPCLSJJ& compare) const {
    bool result= JPCLS::operator<(compare);

    if (result) return result;

    if ( S == compare.S){
      if ( J1 < compare.J1) result=true;
      else if (J1 == compare.J1){
	if ( J2 < compare.J2) result=true;
	else if (J2 == compare.J2){
	  if ( Lambda1 < compare.Lambda1) result=true;
	  else if (Lambda1 == compare.Lambda1){
	    if ( Lambda2 < compare.Lambda2) result=true;
	  }
	}
      }
    }

    return result;
}

  virtual std::string name() const {
    std::string result=JPCLS::name();

    std::stringstream tmpStrStreamJ1;
    tmpStrStreamJ1 << J1;
    std::stringstream tmpStrStreamLam1;
    tmpStrStreamLam1 << Lambda1;
    std::stringstream tmpStrStreamJ2;
    tmpStrStreamJ2 << J2;
    std::stringstream tmpStrStreamLam2;
    tmpStrStreamLam2 << Lambda2;

    result+="J1_"+tmpStrStreamJ1.str()+"Lam1"+tmpStrStreamLam1.str()+"J2_"
      +tmpStrStreamJ2.str()+"Lam2"+tmpStrStreamLam2.str();
    return result;
  }


  virtual void print(std::ostream& os) const {
    JPCLS::print(os);
    os <<"\tJ1,lam1=" << J1 << "," << Lambda1 <<"\tJ2,lam2=" << J2 << "," << Lambda2
       << "\tCGLS=" << CGLS << "\tCGJJ=" << CGJJ << "\tpreFactorAll=" << prefactorAll;
  }
};

struct JPCLSJJls : public JPCLSJJ {
  Spin l;
  Spin s;

  JPCLSJJls(std::shared_ptr<const JPCLSJJ> theJPCLSJJ, const Spin& thel,
	    const Spin& thes): JPCLSJJ(theJPCLSJJ){
    l=thel;
    s=thes;
  }

  JPCLSJJls(std::shared_ptr<const JPCLSJJls> theJPCLSJJls): JPCLSJJ(theJPCLSJJls) {
    l=theJPCLSJJls->l;
    s=theJPCLSJJls->s;
  }

  virtual bool operator==(const jpcRes& compare) const {
    return jpcRes::operator==(compare);
  }

  virtual bool operator==(const JPCLSJJls& compare) const {
    bool result = JPCLSJJ::operator==(compare);
    if ( fabs(l-compare.l)>1e-8 || fabs(s-compare.s)>1e-8 ) result=false;
    return result;
  }

  virtual bool operator<(const JPCLSJJls& compare) const {

    bool result= JPCLSJJ::operator<(compare);

    if (result) return result;

    if ( Lambda2 == compare.Lambda2){
      if ( l < compare.l) result=true;
      else if (l == compare.l){
	if ( s < compare.s) result=true;
      }
    }

    return result;
}

  virtual std::string name() const {
    std::string result=JPCLSJJ::name();

    std::stringstream tmpStrStreaml;
    tmpStrStreaml << l;
    std::stringstream tmpStrStreams;
    tmpStrStreams << s;

    result+="l"+tmpStrStreaml.str()+"s"+tmpStrStreams.str();
    return result;
  }

  virtual void print(std::ostream& os) const {
    JPCLSJJ::print(os);
    os <<"\tl,s=" << l << "," << s;
  }
};


struct JPCLSls : public JPCLS {
  Spin l;
  Spin s;

  JPCLSls(std::shared_ptr<const JPCLS> theJPCLS, const Spin& thel, const Spin& thes):
    JPCLS(theJPCLS) {
    l=thel;
    s=thes;
  }

  JPCLSls(std::shared_ptr<const JPCLSls> theJPCLSls): JPCLS(theJPCLSls) {
    l=theJPCLSls->l;
    s=theJPCLSls->s;
  }

  virtual bool operator==(const jpcRes& compare) const {
    return jpcRes::operator==(compare);
  }

  virtual bool operator==(const JPCLSls& compare) const {
    bool result = JPCLS::operator==(compare);
    if ( fabs(l-compare.l)>1e-8 || fabs(s-compare.s)>1e-8 ) result=false;
    return result;
  }

  virtual bool operator<(const JPCLSls& compare) const {

    bool result= JPCLS::operator<(compare);

    if (result) return result;

    if ( S == compare.S){
      if ( l < compare.l) result=true;
      else if (l == compare.l){
	if ( s < compare.s) result=true;
      }
    }

    return result;
}

  virtual std::string name() const {
    std::string result=JPCLS::name();

    std::stringstream tmpStrStreaml;
    tmpStrStreaml << l;
    std::stringstream tmpStrStreams;
    tmpStrStreams << s;

    result+="l"+tmpStrStreaml.str()+"s"+tmpStrStreams.str();
    return result;
  }


  virtual void print(std::ostream& os) const {
    JPCLS::print(os);
    os <<"\tl,s=" << l << "," << s;
  }
};

struct JPCSM : public jpcRes {
  Spin S;
  Spin M;
  JPCSM(std::shared_ptr<const jpcRes> theJPC, const Spin& theS, const Spin& theM):
    jpcRes(theJPC){
    S=theS;
    M=theM;
  }

  virtual bool operator==(const jpcRes& compare) const {
    return jpcRes::operator==(compare);
  }

  virtual bool operator==(const JPCSM& compare) const {
    bool result=jpcRes::operator==(compare);;
    if ( fabs(S-compare.S)>1e-8 || fabs(M-compare.M)>1e-8 ) result=false;
    return result;
  }

  virtual bool operator<(const JPCSM& compare) const {
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

  virtual std::string name() const {
    std::string result=jpcRes::name();
    std::stringstream tmpStrStreamS;
    tmpStrStreamS << S;
    std::stringstream tmpStrStreamM;
    tmpStrStreamM << M;

    result+="S"+tmpStrStreamS.str()+"M"+tmpStrStreamM.str();
    return result;
  }

  virtual void print(std::ostream& os) const {
    jpcRes::print(os);
    os <<"\tS=" << S << "\tM=" << M
       << std::endl;
  }
};


struct JPCLSM : public JPCLS {
  std::shared_ptr<const jpcRes> jpc;
  Spin M;
  double ClebschG;

  JPCLSM(std::shared_ptr<const JPCLS> theJPCLS, const Spin theM, const double theClebschG) :
    JPCLS(theJPCLS) {
    M=theM;
    ClebschG=theClebschG;
  }

  virtual bool operator==(const jpcRes& compare) const {
    return jpcRes::operator==(compare);
  }

  virtual bool operator==(const JPCLSM& compare) const {
    bool result=jpcRes::operator==(compare);;
    if ( fabs(L-compare.L)>1e-8 || fabs(S-compare.S)>1e-8 || fabs(M-compare.M)>1e-8 )
      result=false;
    return result;
  }

  virtual bool operator<(const JPCLSM& compare) const {
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


  virtual std::string name() const {
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


  virtual void print(std::ostream& os) const {
    JPCLS::print(os);
    os <<"\tM=" << M
       <<"\tClebschGordan=" << ClebschG
       << std::endl;
  }
};


template<typename T>
void fillVec(std::shared_ptr<const T> currentRes,
	     std::vector< std::shared_ptr<const T> >& theVec) {

  const T* current=currentRes.get();
  if (0==current) {
      std::cout << "ALERT: shared object contains 0 pointer!!!" << std::endl;
    exit(1);
  }

  typename std::vector< std::shared_ptr<const T> >::const_iterator it;
  bool found=false;
  for ( it=theVec.begin(); it!=theVec.end(); ++it){
    const T* resIt=(*it).get();
    if (0==resIt) {
      std::cout << "ALERT: shared object contains 0 pointer!!!" << std::endl;
      exit(1);
    }
    if ( (*resIt) == (*current) ){
      found=true;
      continue;
    }
  }
  if (!found) theVec.push_back(currentRes);
}


template<typename T1, typename T2>
std::vector< std::shared_ptr<const T2 > > extractStates (std::vector< std::shared_ptr<const T1 > >& vec1, std::vector< std::shared_ptr<const T2 > >& vec2) {

  typename std::vector< std::shared_ptr<const T2 > > result;
  typename std::vector< std::shared_ptr<const T1 > >::const_iterator it1;
  typename std::vector< std::shared_ptr<const T2 > >::const_iterator it2;

  for ( it1=vec1.begin(); it1!=vec1.end(); ++it1) {
    for ( it2 = vec2.begin(); it2!= vec2.end(); ++it2) {
      if ( *(*it1)== *(*it2) )  fillVec( (*it2), result);
    }
  }

   return result;
}

template<typename T1, typename T2, typename T3>
void fillStatesInitialDec( std::vector< std::shared_ptr<const T1 > >& initial,
			   std::vector< std::shared_ptr<const T2 > >& dec,
			   std::vector< std::shared_ptr<const T3 > >& initialDecToFill) {

  typename std::vector< std::shared_ptr<const T1 > >::const_iterator it1;
  typename std::vector< std::shared_ptr<const T2 > >::const_iterator it2;

  for ( it1=initial.begin(); it1!=initial.end(); ++it1) {

    jpcRes thejpcRes( (*it1)->J, (*it1)->P, (*it1)->C );
    thejpcRes.print(std::cout);
    for ( it2=dec.begin(); it2!=dec.end(); ++it2){
      if( *(*it2) == thejpcRes){
	std::shared_ptr<const T3> currentInitialDecState(new T3( (*it1), (*it2)->L,
								   (*it2)->S));
	initialDecToFill.push_back(currentInitialDecState);
      }
    }
  }
}


void validJPCLS(std::shared_ptr<const jpcRes> motherRes, std::shared_ptr<const jpcRes> daughterRes1, std::shared_ptr<const jpcRes> daughterRes2, std::vector< std::shared_ptr<const JPCLS> >& theJPCLSVec);

class Particle;

void validJPCLS(std::shared_ptr<const jpcRes> motherRes, Particle* daughter1, Particle* daughter2, std::vector< std::shared_ptr<const JPCLS> >& theJPCLSVec,  bool useCParity=true, int gParityMother=0, bool useIsospin=false);

void validJPClamlam(std::shared_ptr<const jpcRes> motherRes, Particle* daughter1, Particle* daughter2, std::vector< std::shared_ptr<const JPClamlam> >& theJPClamlamVec,  bool useCParity=true, int gParityMother=0, bool useIsospin=false);

std::shared_ptr<jpcRes> getJPCPtr(Particle* theParticle);
std::shared_ptr<IGJPC> getIGJPCPtr(Particle* theParticle);

