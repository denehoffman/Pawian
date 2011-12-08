#include "qft++/matrix/KpoleMatrix.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"

KpoleMatrix::KpoleMatrix(const vector<double>& g_i, const double mass_0, const vector< pair<double, double> >& mDecPair, int orbMom):
  Matrix< complex<double> >::Matrix(int(g_i.size()), int(g_i.size()))
  , _orbMom(orbMom)
  , _g_i(g_i)
  , _poleMass(mass_0)
  ,_mDecPair(mDecPair) 
{
  //calculate phase-space factors

  for (int i=0; i<int(_g_i.size()); ++i){
    complex<double> tmpRhoPoleFactors=phaseSpaceFac(_poleMass,_mDecPair[i].first, _mDecPair[i].second); 
    _rhoPoleFactors.push_back(tmpRhoPoleFactors);
    
    complex<double> tmpBreakupFactors=breakupMomQ(_poleMass,_mDecPair[i].first, _mDecPair[i].second);
    _breakupFactors.push_back(tmpBreakupFactors);
    
  }  
}

KpoleMatrix::KpoleMatrix(const KpoleMatrix &theCopy):
  Matrix< complex<double> >::Matrix(int(theCopy._g_i.size()), int(theCopy._g_i.size()))
  , _orbMom(theCopy._orbMom)
  ,_g_i(theCopy._g_i)
  , _poleMass(theCopy._poleMass)
   ,_mDecPair(theCopy._mDecPair) 
{
  //calculate phase-space factors

  for (int i=0; i<int(_g_i.size()); ++i){
    complex<double> tmpRhoPoleFactors=phaseSpaceFac(_poleMass,_mDecPair[i].first, _mDecPair[i].second); 
    _rhoPoleFactors.push_back(tmpRhoPoleFactors);
    
    complex<double> tmpBreakupFactors=breakupMomQ(_poleMass,_mDecPair[i].first, _mDecPair[i].second);
    _breakupFactors.push_back(tmpBreakupFactors);
    
  }
}

KpoleMatrix::~KpoleMatrix(){

}

void KpoleMatrix::updateMatrix(const double mass){
  //calculate phase-space factors


  vector< complex<double> > rhoFactors;
  vector< complex<double> > barrierFactors;

  for (int i=0; i< int(_g_i.size()); ++i){
    complex<double> tmpRhoFactors=phaseSpaceFac(mass,_mDecPair[i].first, _mDecPair[i].second);
    rhoFactors.push_back(tmpRhoFactors);

    complex<double> tmpBarrierFactor(1.,0.);    
    if(_orbMom>0) tmpBarrierFactor=pow(breakupMomQ(mass, _mDecPair[i].first, _mDecPair[i].second)/_breakupFactors[i], _orbMom);
    barrierFactors.push_back(tmpBarrierFactor);
  }


  for (int i=0; i< int(_g_i.size()); ++i){
    for (int j=0; j< int(_g_i.size()); ++j){
      this->operator()(i,j)= ( _g_i[i]*sqrt(conj(rhoFactors[i]))*barrierFactors[i]*_g_i[j]*sqrt(rhoFactors[j])*barrierFactors[j])/(_poleMass*_poleMass-mass*mass);
 //        this->operator()(i,j)= ( _g_i[i]*sqrt(rhoFactors[i]*rhoFactors[j])*_g_i[j])/(_poleMass*_poleMass-mass*mass);
//       this->operator()(i,j)= (( _g_i[i]*_g_i[j])/(_poleMass*_poleMass-mass*mass))*sqrt(rhoFactors[i]*rhoFactors[j]);
    }
  }
  
}

void KpoleMatrix::updateMatrixRel(const double mass){

  updateMatrix(mass);

  vector< complex<double> > rhoFactors;
  for (int i=0; i< int(_g_i.size()); ++i){
    complex<double> tmpRhoFactors=phaseSpaceFac(mass,_mDecPair[i].first, _mDecPair[i].second);
    rhoFactors.push_back(tmpRhoFactors);
  }

  for (int i=0; i< int(_g_i.size()); ++i){
    for (int j=0; j< int(_g_i.size()); ++j){
//       this->operator()(i,j) /= sqrt(_rhoPoleFactors[i]*_rhoPoleFactors[j]);
      this->operator()(i,j) /= sqrt(conj(rhoFactors[i])*rhoFactors[j]);
     }
   }
  
}
