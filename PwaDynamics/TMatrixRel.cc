#include "PwaDynamics/TMatrixRel.hh"
#include "PwaDynamics/KMatrixBase.hh"
#include "qft++/matrix/IdentityMatrix.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"


TMatrixRel::TMatrixRel(boost::shared_ptr<KMatrixBase> Kmatrix) :
  TMatrixBase(Kmatrix)
 {
 }


TMatrixRel::~TMatrixRel(){
}

void TMatrixRel::evalMatrix(const double mass){
  TMatrixBase::evalMatrix(mass);

  vector<boost::shared_ptr<AbsPhaseSpace> > phpVec=_Kmatrix->phaseSpaceVec();
  Matrix< complex<double> > theRhoMatrix(NumRows(),NumRows());

  for (int i=0; i<NumRows(); ++i){
    for (int j=0; j<NumRows(); ++j){
      theRhoMatrix(i,j)=complex<double> (0.,0.);
      if (i==j){
        theRhoMatrix(i,j) = phpVec[j]->factor(mass); 
      }
    }
  }

  complex<double> imagCompl(0.,1.);  
  IdentityMatrix< complex<double> > theIdMatrix(NumRows());
  Matrix< complex< double > > tmpDenomMatrCompl = theIdMatrix-imagCompl*(*_Kmatrix)*theRhoMatrix;
  
  Matrix< complex< double > > tmpDenomMatrInv=tmpDenomMatrCompl;
  tmpDenomMatrInv.invert();
  
  Matrix< complex <double> > currentTMatr=tmpDenomMatrInv*(*_Kmatrix);
  
  for (int i=0; i<currentTMatr.NumRows(); ++i){
    for (int j=0; j<currentTMatr.NumCols(); ++j){
      this->operator()(i,j)=currentTMatr(i,j);
    }
  }
}
