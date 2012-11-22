#include "PwaDynamics/PVectorSlowCorRel.hh"
#include "PwaDynamics/PPole.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++/matrix/IdentityMatrix.hh"
#include "ErrLogger/ErrLogger.hh"

PVectorSlowCorRel::PVectorSlowCorRel(vector<boost::shared_ptr<PPole> > Ppoles, vector<boost::shared_ptr<AbsPhaseSpace> > phpVecs, std::vector< complex<double> >& fProdVec, double s0prod) :
  PVectorRel(Ppoles, phpVecs)
  ,_fProdVec(fProdVec)
  ,_s0prod(s0prod)
 {
 } 

PVectorSlowCorRel::~PVectorSlowCorRel(){
}

void PVectorSlowCorRel::evalMatrix(const double mass){

  Matrix< complex<double> > thePVector(NumRows(), 1);
  vector<boost::shared_ptr<PPole> >::iterator it;
  for (it =_Ppoles.begin(); it != _Ppoles.end(); ++it){
    (*it)->evalMatrix(mass);
    thePVector += *(*it);
  }

  for (int i=0; i<thePVector.NumRows(); ++i){  
    this->operator()(i,0)=thePVector(i,0)+ _fProdVec[i]*(1.0-_s0prod)/(mass*mass-_s0prod);
  } 

}

void PVectorSlowCorRel::updateFprod(int i, complex<double> fProd) {
  if ( i >= int(_fProdVec.size()) ){
    Alert << "no fProd parameter with index " << i << " available!!!" << endmsg;
    exit(0);
  }
  _fProdVec[i]=fProd;
}
