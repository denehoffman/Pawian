#include "PwaDynamics/FVector.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"


FVector::FVector(boost::shared_ptr<KMatrixBase> Kmatrix, boost::shared_ptr<PVectorRel> Pvector) :
  Matrix< complex<double> >::Matrix(Pvector->NumRows(), Pvector->NumCols())
  , _Kmatrix(Kmatrix)
  ,_Pvector(Pvector)
 {
 }

FVector::FVector(int numRows) :
  Matrix< complex<double> >::Matrix(numRows, 1)
 {
 }

FVector::~FVector(){
}

void FVector::evalMatrix(const double mass){
  _Kmatrix->evalMatrix(mass);
  _Pvector->evalMatrix(mass);

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
  
  Matrix< complex <double> > currentTMatr=tmpDenomMatrInv*(*_Pvector);

  for (int i=0; i<currentTMatr.NumRows(); ++i){
    this->operator()(i,0)=currentTMatr(i,0);
  }

}

