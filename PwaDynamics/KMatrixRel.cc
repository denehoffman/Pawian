#include "PwaDynamics/KMatrixRel.hh"
#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++/matrix/IdentityMatrix.hh"

KMatrixRel::KMatrixRel(vector<boost::shared_ptr<KPole> > Kpoles, vector<boost::shared_ptr<AbsPhaseSpace> > phpVecs) :
  KMatrixBase(Kpoles, phpVecs)
 {
 }


KMatrixRel::~KMatrixRel(){
}

void KMatrixRel::evalMatrix(const double mass){

  Matrix< complex<double> > theKMatrix(NumRows(), NumRows());
  vector<boost::shared_ptr<KPole> >::iterator it;
  for (it =_KPoles.begin(); it != _KPoles.end(); ++it){
    (*it)->evalMatrix(mass);
    theKMatrix += *(*it);
  }

  for (int i=0; i<theKMatrix.NumRows(); ++i){
    for (int j=0; j<theKMatrix.NumCols(); ++j){
      this->operator()(i,j)=theKMatrix(i,j);
    }
  }
}
