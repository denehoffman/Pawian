#ifndef _FitParamErrorMatrixStreamer_H
#define _FitParamErrorMatrixStreamer_H

#include <string>
#include <vector>


class FitParamErrorMatrixStreamer {
  
public :
  FitParamErrorMatrixStreamer(  std::string errFile  );
  virtual ~FitParamErrorMatrixStreamer(){};
  
  void matrixData( std::vector<double>& theData, int& nrows  ){ 
    theData=_matrixData;  
    int sizeOfData=(int) theData.size();
    nrows=1;
    while(nrows*nrows!=sizeOfData){
      nrows++;
    
    }
  }
  
private :
  std::vector<double> _matrixData;
};

#endif

