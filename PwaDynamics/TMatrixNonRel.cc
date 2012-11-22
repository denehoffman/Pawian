#include "PwaDynamics/TMatrixNonRel.hh"
#include "PwaDynamics/KMatrixNonRel.hh"
#include "qft++/matrix/IdentityMatrix.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"


TMatrixNonRel::TMatrixNonRel(boost::shared_ptr<KMatrixNonRel> Kmatrix) :
  TMatrixBase(Kmatrix)
 {
 }


TMatrixNonRel::~TMatrixNonRel(){
}

void TMatrixNonRel::evalMatrix(const double mass){
  TMatrixBase::evalMatrix(mass);

  complex<double> imagCompl(0.,1.);  
  IdentityMatrix< complex<double> > theIdMatrix(NumRows());
  Matrix< complex< double > > tmpDenomMatrCompl = theIdMatrix-imagCompl*(*_Kmatrix);
  
  Matrix< complex< double > > tmpDenomMatrInv=tmpDenomMatrCompl;
  tmpDenomMatrInv.invert();
  
  Matrix< complex <double> > currentTMatr=tmpDenomMatrInv*(*_Kmatrix);
  
  for (int i=0; i<currentTMatr.NumRows(); ++i){
    for (int j=0; j<currentTMatr.NumCols(); ++j){
      this->operator()(i,j)=currentTMatr(i,j);
    }
  }
}
