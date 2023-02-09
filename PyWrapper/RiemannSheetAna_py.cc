#include "PyWrapper/RiemannSheetAna_py.hh"
#include "KMatrixExtract/TMatrixExtrBase.hh"
#include "KMatrixExtract/TMatrixResidueExtr.hh"
#include "KMatrixExtract/FVectorResidueExtr.hh"
#include "ConfigParser/pipiScatteringParser.hh"

RiemannSheetAna_py::RiemannSheetAna_py()
{
  char* argv[7];
  argv[0]=(char*)"$TOP_DIR/bin/TMatrixGeneralApp";
  argv[1]=(char*)"-c";
  argv[2]=(char*)"/data/duldul/bertram/PawianPyWrapper/Pawian/Examples/KMatrixExtract/scatteringa2PhasePiEta.cfg";
  argv[3]=(char*)"--paramFile";
  argv[4]=(char*)"/data/duldul/bertram/PawianPyWrapper/Pawian/Examples/KMatrixExtract/benchmarkResult.dat";
  argv[5]=(char*)"--serializationFile";
  argv[6]=(char*)"/data/duldul/bertram/PawianPyWrapper/Pawian/Examples/KMatrixExtract/serializedOutput.dat";
  _pipiParser = new pipiScatteringParser(7, argv);
  _tMatrixExtrBasePtr = new TMatrixExtrBase(_pipiParser);
  _tMatrixResidueExtrPtr = new TMatrixResidueExtr(_pipiParser);

  argv[2]=(char*)"/data/duldul/bertram/PawianPyWrapper/Pawian/Examples/KMatrixExtract/scatteringa2FVectorIntensityPiEta.cfg";

  _pipiParserFVec = new pipiScatteringParser(7, argv);
  _fVectorResidueExtrPtr = new FVectorResidueExtr(_pipiParserFVec);
}

RiemannSheetAna_py::~RiemannSheetAna_py(){
 delete _pipiParser;
 delete _pipiParserFVec;
 delete _tMatrixExtrBasePtr;
 delete _tMatrixResidueExtrPtr;
 delete _fVectorResidueExtrPtr;
}

double RiemannSheetAna_py::calcTMat(double eReal, double eImag){
  return _tMatrixExtrBasePtr->calcTMatrix(eReal, eImag);
}

void RiemannSheetAna_py::calcResidue(){
  _tMatrixResidueExtrPtr->Calculation();
}

void RiemannSheetAna_py::calcFVecResidue(){
  _fVectorResidueExtrPtr->Calculation();
}

void RiemannSheetAna_py::SetParamValue(const std::string & paramName, double paramVal){
  _tMatrixExtrBasePtr->SetParamValue(paramName, paramVal);
}



