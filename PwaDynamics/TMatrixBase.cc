#include "PwaDynamics/TMatrixBase.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"


TMatrixBase::TMatrixBase(boost::shared_ptr<KMatrixBase> Kmatrix) :
  Matrix< complex<double> >::Matrix(int(Kmatrix->phaseSpaceVec().size()), int(Kmatrix->phaseSpaceVec().size()))
  , _Kmatrix(Kmatrix)
 {
 }


TMatrixBase::~TMatrixBase(){
}

void TMatrixBase::evalMatrix(const double mass){
  _Kmatrix->evalMatrix(mass);
}

