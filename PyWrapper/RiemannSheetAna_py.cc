#include "PyWrapper/RiemannSheetAna_py.hh"
#include "KMatrixExtract/TMatrixExtrBase.hh"
#include "ConfigParser/pipiScatteringParser.hh"

RiemannSheetAna_py::RiemannSheetAna_py()
{
  char* argv[5];
  argv[0]=(char*)"$TOP_DIR/bin/TMatrixGeneralApp";
  argv[1]=(char*)"-c";
  argv[2]=(char*)"/data/duldul/bertram/PawianPyWrapper/Pawian/Examples/KMatrixExtract/scatteringa2PhasePiEta.cfg";
  argv[3]=(char*)"--paramFile";
  argv[4]=(char*)"/data/duldul/bertram/PawianPyWrapper/Pawian/Examples/KMatrixExtract/benchmarkResult.dat";
  _pipiParser = new pipiScatteringParser(5, argv);
  _tMatrixExtrBasePtr = new TMatrixExtrBase(_pipiParser);
}

RiemannSheetAna_py::~RiemannSheetAna_py(){
 delete _pipiParser;
 delete _tMatrixExtrBasePtr;
}

double RiemannSheetAna_py::calcTMat(double eReal, double eImag){
  return _tMatrixExtrBasePtr->calcTMatrix(eReal, eImag);
}




