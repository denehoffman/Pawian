#include "PwaDynamics/KMatrixBase.hh"
#include "PwaDynamics/KPole.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++/matrix/IdentityMatrix.hh"

KMatrixBase::KMatrixBase(vector<boost::shared_ptr<KPole> > Kpoles, vector<boost::shared_ptr<AbsPhaseSpace> > phpVecs) :
  Matrix< complex<double> >::Matrix(int(phpVecs.size()), int(phpVecs.size()))
  ,_KPoles(Kpoles)
  ,_phpVecs(phpVecs)
 {
 }


KMatrixBase::~KMatrixBase(){
}

