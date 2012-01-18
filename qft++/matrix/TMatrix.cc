#include "qft++/matrix/TMatrix.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++/matrix/IdentityMatrix.hh"

TMatrix::TMatrix(const vector<KpoleMatrix>& theKpoles):
  Matrix< complex<double> >::Matrix(theKpoles[0].NumRows(), theKpoles[0].NumCols())
  ,_Kpoles(theKpoles)
 {
   _currentRhoFactors.resize(theKpoles[0].NumRows());
 }


TMatrix::~TMatrix(){

}

void TMatrix::updateMatrix(const double mass){

  Matrix< complex<double> > theKMatrix(NumRows(), NumRows());
  vector<KpoleMatrix>::iterator it;

  for (it =_Kpoles.begin(); it != _Kpoles.end(); ++it){
    it->updateMatrix(mass);
    theKMatrix += (*it);
  }

  cout << "theKMatrix: " << theKMatrix << endl;
  
  
  complex<double> imagCompl(0.,1.);
  IdentityMatrix< complex<double> > theIdMatrix(NumRows());

  Matrix< complex< double > > tmpDenomMatrCompl = theIdMatrix - imagCompl*theKMatrix;
  
  
  Matrix< complex< double > > tmpDenomMatrInv=tmpDenomMatrCompl;
  tmpDenomMatrInv.invert();
  
  
  //  Matrix< complex< double > > testMatr=tmpDenomMatrCompl*tmpDenomMatrInv;
//   cout << "testMatr: " << testMatr << endl;
  
  
  Matrix< complex <double> > currentTMatr=tmpDenomMatrInv*theKMatrix;
  
  for (int i=0; i<currentTMatr.NumRows(); ++i){
    for (int j=0; j<currentTMatr.NumCols(); ++j){
      
      this->operator()(i,j)=currentTMatr(i,j);
      
    }
  }
  
  
}

void TMatrix::updateMatrixRel(const double mass){

  Matrix< complex<double> > theKMatrix(NumRows(), NumRows());
  vector<KpoleMatrix>::iterator it;

  for (it =_Kpoles.begin(); it != _Kpoles.end(); ++it){
     it->updateMatrixRel(mass);
     theKMatrix += (*it);
  }

  //  cout << "theKMatrixRel: " << theKMatrix << endl;

  Matrix< complex<double> > theRhoMatrix(NumRows(),NumRows());

  for (int i=0; i<NumRows(); ++i){
    for (int j=0; j<NumRows(); ++j){
      theRhoMatrix(i,j)=complex<double> (0.,0.);
      if (i==j){
	//          complex<double> tmpRhoPoleFactors=phaseSpaceFac(mass,_mDecPair[j].first, _mDecPair[j].second);
	// theRhoMatrix(i,j) = tmpRhoPoleFactors;
	complex<double> tmpRhoPoleFactors=phaseSpaceFac(mass,_Kpoles[0].decayPair(j).first, _Kpoles[0].decayPair(j).second);
	_currentRhoFactors[i]=tmpRhoPoleFactors;
	theRhoMatrix(i,j) = tmpRhoPoleFactors; 
       }
    }
  }


  complex<double> imagCompl(0.,1.);
  IdentityMatrix< complex<double> > theIdMatrix(NumRows());

  Matrix< complex< double > > tmpDenomMatrCompl = theIdMatrix - imagCompl*theKMatrix*theRhoMatrix;

  Matrix< complex< double > > tmpDenomMatrInv=tmpDenomMatrCompl;
  tmpDenomMatrInv.invert();


  //  Matrix< complex< double > > testMatr=tmpDenomMatrCompl*tmpDenomMatrInv;
  //   cout << "testMatrRel: " << testMatr << endl;



   Matrix< complex <double> > currentTMatr=tmpDenomMatrInv*theKMatrix;


  for (int i=0; i<currentTMatr.NumRows(); ++i){
    for (int j=0; j<currentTMatr.NumCols(); ++j){
      
      this->operator()(i,j)=currentTMatr(i,j);
    }
  }


//   if (NumRows()==2){
//     Matrix< complex< double > > theTestTMat(NumRows(),NumRows());

//      complex<double> theTestD=_Kpoles[0](0,0)*_Kpoles[0](1,1)-_Kpoles[0](0,1)*_Kpoles[0](0,1);
 
// //     complex<double> theTestD(0.,0.);
// //     cout << "theTestDRel: " << theTestD << endl;
//     complex<double> theTestDenom=1.-theRhoMatrix(0,0)*theRhoMatrix(1,1)*theTestD - imagCompl*(theRhoMatrix(0,0)*_Kpoles[0](0,0)+theRhoMatrix(1,1)*_Kpoles[0](1,1));


//     theTestTMat(0,0)=(_Kpoles[0](0,0) - imagCompl*theRhoMatrix(1,1)*theTestD)/theTestDenom;
//     theTestTMat(0,1)=_Kpoles[0](0,1)/theTestDenom;
//     theTestTMat(1,0)=_Kpoles[0](1,0)/theTestDenom;
//     theTestTMat(1,1)=(_Kpoles[0](1,1) - imagCompl*theRhoMatrix(0,0)*theTestD)/theTestDenom;

//  //    cout << "theTestTMatRel: " << theTestTMat << endl;

//   }

}
