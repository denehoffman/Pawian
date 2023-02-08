#include <iostream>
#include <math.h>
#include <boost/python.hpp>

using namespace boost::python;
class TMatrixExtrBase;
class pipiScatteringParser;

class RiemannSheetAna_py {

public:
  
  RiemannSheetAna_py();
  ~RiemannSheetAna_py();
  double calcTMat(double eReal, double eImag);

private:
  pipiScatteringParser* _pipiParser;
  TMatrixExtrBase* _tMatrixExtrBasePtr;
};

BOOST_PYTHON_MODULE(RiemannSheetAna_py)
{
    class_<RiemannSheetAna_py>("RiemannSheetAna_py")
      .def("calcTMat", &RiemannSheetAna_py::calcTMat)
      ;
}

