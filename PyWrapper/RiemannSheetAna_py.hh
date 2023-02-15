#include <iostream>
#include <math.h>
#include <boost/python.hpp>

using namespace boost::python;

class TMatrixExtrBase;
class pipiScatteringParser;
class TMatrixResidueExtr;
class FVectorResidueExtr;

class RiemannSheetAna_py {

public:
  
  RiemannSheetAna_py();
  RiemannSheetAna_py(boost::python::list);
  ~RiemannSheetAna_py();
  double calcTMat(double eReal, double eImag);
  void calcResidue();
  void calcFVecResidue();
  void SetParamValue(const std::string & paramName, double paramVal);
  void SetSheet(const std::string & sheetName);
  
private:
  pipiScatteringParser* _pipiParser;
  pipiScatteringParser* _pipiParserFVec;
  TMatrixExtrBase* _tMatrixExtrBasePtr;
  TMatrixResidueExtr* _tMatrixResidueExtrPtr;
  FVectorResidueExtr* _fVectorResidueExtrPtr;
};

BOOST_PYTHON_MODULE(RiemannSheetAna_py)
{
    class_<RiemannSheetAna_py>("RiemannSheetAna_py")
      .def(init<boost::python::list>())
      .def("calcTMat", &RiemannSheetAna_py::calcTMat)
      .def("calcResidue", &RiemannSheetAna_py::calcResidue)
      .def("calcFVecResidue", &RiemannSheetAna_py::calcFVecResidue)
      .def("SetParamValue", &RiemannSheetAna_py::SetParamValue)
      .def("SetSheet", &RiemannSheetAna_py::SetSheet)
      ;
}

