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
  void SetParamValue(const std::string & paramName, double paramVal);
  double GetParamValue(const std::string & paramName);
  void SetSheet(const std::string & sheetName);
  int sheetSize();
  
private:
  pipiScatteringParser* _pipiParser;
  TMatrixExtrBase* _tMatrixExtrBasePtr;
  TMatrixResidueExtr* _tMatrixResidueExtrPtr;
};

BOOST_PYTHON_MODULE(RiemannSheetAna_py)
{
  class_<RiemannSheetAna_py>("RiemannSheetAna_py")
    .def(init<>())
    .def(init<boost::python::list>())
    .def("calcTMat", &RiemannSheetAna_py::calcTMat)
    .def("calcResidue", &RiemannSheetAna_py::calcResidue)
    .def("SetParamValue", &RiemannSheetAna_py::SetParamValue)
    .def("GetParamValue", &RiemannSheetAna_py::GetParamValue)
    .def("SetSheet", &RiemannSheetAna_py::SetSheet)
    .def("sheetSize", &RiemannSheetAna_py::sheetSize)
    ;
}

