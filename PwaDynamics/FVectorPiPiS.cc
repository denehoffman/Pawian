#include "PwaDynamics/FVectorPiPiS.hh"
#include "PwaDynamics/KMatrixPiPiS.hh"
#include "PwaDynamics/PVectorRel.hh"
#include "PwaDynamics/PPole.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include <boost/multi_array.hpp>

FVectorPiPiS::FVectorPiPiS() :
  FVector(5)
 {
   boost::shared_ptr<KMatrixPiPiS> theKMatrix(new KMatrixPiPiS());
   vector<boost::shared_ptr<AbsPhaseSpace> > phpVecs=theKMatrix->phaseSpaceVec();
   vector<boost::shared_ptr<KPole> > kPoles=theKMatrix->kpoles();

   vector<boost::shared_ptr<PPole> > thePpoles;
   complex<double> defaultBeta(1.,0.); 
     
   vector<boost::shared_ptr<KPole> >::iterator it;
   for (it=kPoles.begin(); it!=kPoles.end(); ++it){ 
     std::vector<double> currentGFactors=(*it)->gFactors();
     boost::shared_ptr<PPole> currentPPole(new PPole(defaultBeta, currentGFactors, (*it)->poleMass()));
     thePpoles.push_back(currentPPole);     
   } 

   std::vector<complex <double> > fProdVec;

   for (int i=0; i<int(phpVecs.size()); ++i){
     complex<double> currentVal(1.0,0.);
     fProdVec.push_back(currentVal);
  }

  double s0Prod=-0.0737;
  boost::shared_ptr<PVectorSlowCorRel> thePVector(new PVectorSlowCorRel(thePpoles, phpVecs, fProdVec, s0Prod));   
  _Kmatrix=theKMatrix;
  _Pvector=thePVector;
  _pVectorCor=thePVector;
     
 }


FVectorPiPiS::~FVectorPiPiS(){
}

