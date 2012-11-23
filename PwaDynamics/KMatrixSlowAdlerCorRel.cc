#include "PwaDynamics/KMatrixSlowAdlerCorRel.hh"
#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++/matrix/IdentityMatrix.hh"

KMatrixSlowAdlerCorRel::KMatrixSlowAdlerCorRel(vector<boost::shared_ptr<KPole> > Kpoles, vector<boost::shared_ptr<AbsPhaseSpace> > phpVecs, boost::shared_ptr<array_type_2d> fscatPtr, double s0_scat, double sAdler0, double sAdler) :
  KMatrixBase(Kpoles, phpVecs)
  ,_fScatPtr(fscatPtr)
  ,_s0Scat(s0_scat)
  ,_sAdler0(sAdler0)
  ,_sAdler(sAdler)
 {
 }

KMatrixSlowAdlerCorRel::KMatrixSlowAdlerCorRel(int dim) :
  KMatrixBase(dim, dim)
  ,_s0Scat(0.)
  ,_sAdler0(0.)
  ,_sAdler(0.)
{
}

KMatrixSlowAdlerCorRel::~KMatrixSlowAdlerCorRel(){
}

void KMatrixSlowAdlerCorRel::evalMatrix(const double mass){

  Matrix< complex<double> > theKMatrix(NumRows(), NumRows());
  vector<boost::shared_ptr<KPole> >::iterator it;
  for (it =_KPoles.begin(); it != _KPoles.end(); ++it){
    (*it)->evalMatrix(mass);
    theKMatrix += *(*it);
  }

  for (int i=0; i<theKMatrix.NumRows(); ++i){
    for (int j=0; j<theKMatrix.NumCols(); ++j){
      this->operator()(i,j)=theKMatrix(i,j)+complex<double> ((*_fScatPtr)[i][j]*(1.0-_s0Scat)/(mass*mass-_s0Scat), 0.);
      this->operator()(i,j)*=(1.-_sAdler0)/(mass*mass-_sAdler0)*(mass*mass-_sAdler*0.1349766*0.1349766/2.);
    }
  }
}
