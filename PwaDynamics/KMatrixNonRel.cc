#include "PwaDynamics/KMatrixNonRel.hh"
#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++/matrix/IdentityMatrix.hh"

KMatrixNonRel::KMatrixNonRel(vector<boost::shared_ptr<KPole> > Kpoles, vector<boost::shared_ptr<AbsPhaseSpace> > phpVecs) :
  KMatrixBase(Kpoles, phpVecs)
 {
 }


KMatrixNonRel::~KMatrixNonRel(){
}

void KMatrixNonRel::evalMatrix(const double mass){

  vector< complex<double> > rhoFactors;

  for (int i=0; i< int(_phpVecs.size()); ++i){
    rhoFactors.push_back(_phpVecs[i]->factor(mass));
  }

  Matrix< complex<double> > theKMatrix(NumRows(), NumRows());
  vector<boost::shared_ptr<KPole> >::iterator it;
  for (it =_KPoles.begin(); it != _KPoles.end(); ++it){
    (*it)->evalMatrix(mass);
    theKMatrix += *(*it);
  }

  for (int i=0; i<theKMatrix.NumRows(); ++i){
    for (int j=0; j<theKMatrix.NumCols(); ++j){
      this->operator()(i,j)=theKMatrix(i,j)*sqrt(rhoFactors[i]*rhoFactors[j]);
    }
  }
}
