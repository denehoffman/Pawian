#include "Examples/JpsiGamKsKlKK/FitParamErrorMatrix.hh"

using namespace ROOT::Minuit2;

FitParamErrorMatrix::FitParamErrorMatrix( ROOT::Minuit2::MnUserCovariance covMatrix,  ROOT::Minuit2::MnUserParameters usrParams){
  
  MnUserTransformation trafo = usrParams.Trafo(); 
  _size = usrParams.Parameters().size();
  
  _errMatrix = boost_matrix(_size,_size);
  
  for(unsigned int irow=0; irow<_size; irow++){
    for(unsigned int icol=0; icol<_size; icol++){
      _errMatrix(irow,icol)=0.0;
    }
  }
  
  for(unsigned int row=0;row<covMatrix.Nrow(); row++){
    for(unsigned int col=0;col<covMatrix.Nrow(); col++){
      double element = covMatrix(row,col);
      int trow = trafo.ExtOfInt(row);
      int tcol = trafo.ExtOfInt(col);
      _errMatrix(trow,tcol)=element;
    }
  }
}






FitParamErrorMatrix::~FitParamErrorMatrix(  ){
}


FitParamErrorMatrix::FitParamErrorMatrix(std::vector<double> &theData, int size){
  _size=size;
  _errMatrix = boost_matrix(_size,_size);
  
  unsigned int irow(0), icol(0);
  std::vector<double>::iterator dataIter;
  for(dataIter=theData.begin(); dataIter!=theData.end(); dataIter++){
    double val = *dataIter;
    _errMatrix(irow,icol) = val;
    
    if(icol==_size-1){
      icol=0;
      irow++;
    }else{
      icol++;
    }
  }
}




void FitParamErrorMatrix::Print(std::ostream &os){
  for(unsigned int irow=0;irow<_size;irow++){
    for(unsigned int icol=0;icol<_size;icol++){
      os << _errMatrix(irow,icol) << " ";
    }
    os << std::endl;
  }
}



std::vector<double> FitParamErrorMatrix::Data(){
  std::vector<double> theData;
  for(unsigned int irow=0;irow<_size;irow++){
      for(unsigned int icol=0;icol<_size;icol++){
	theData.push_back(_errMatrix(irow,icol) );
      }
  }
  return theData;


  

}

void FitParamErrorMatrix::Write(std::ostream& os){
  std::vector<double> data = Data();
  std::vector<double>::iterator dataIter = data.begin();
  for(; dataIter!=data.end();dataIter++){
    os << *dataIter << std::endl;
  }
}

