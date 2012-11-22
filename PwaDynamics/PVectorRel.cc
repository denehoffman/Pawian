#include "PwaDynamics/PVectorRel.hh"
#include "PwaDynamics/PPole.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++/matrix/IdentityMatrix.hh"
#include "ErrLogger/ErrLogger.hh"

PVectorRel::PVectorRel(vector<boost::shared_ptr<PPole> > Ppoles, vector<boost::shared_ptr<AbsPhaseSpace> > phpVecs) :
  KMatrixBase(phpVecs, int(phpVecs.size()), 1)
  ,_Ppoles(Ppoles)
 {
 }


PVectorRel::~PVectorRel(){
}

void PVectorRel::evalMatrix(const double mass){

  Matrix< complex<double> > thePVector(NumRows(), 1);
  vector<boost::shared_ptr<PPole> >::iterator it;
  for (it =_Ppoles.begin(); it != _Ppoles.end(); ++it){
    (*it)->evalMatrix(mass);
    thePVector += *(*it);
  }

  for (int i=0; i<thePVector.NumRows(); ++i){
      this->operator()(i,0)=thePVector(i,0);
  }



}

void PVectorRel::updateBeta(int i, complex<double> beta){
  if ( i<0 || i>=int(_Ppoles.size()) ){
    Alert << "no pole with index " << i << " available!!!" << endmsg;
    exit(0);
  }
  _Ppoles[i]->updateBeta(beta);
}
