#include <iostream>
#include <math.h>
#include <boost/python.hpp>
#include "qft++/tensor/Vector4.hh"
//class Vector4;

using namespace std;
using namespace boost::python;

class HeliVec_py {

public:
  
  HeliVec_py();
  ~HeliVec_py(){};
  void calcHeliVec();

private:
  Vector4<double> _motherRefQftVec4;
  Vector4<double> _refQftVec4;
  Vector4<double> _motherQftVec4;
  Vector4<double> _daughterQftVec4;
};

BOOST_PYTHON_MODULE(HeliVec_py)
{
    class_<HeliVec_py>("HeliVec_py")
        .def("calcHeliVec", &HeliVec_py::calcHeliVec)
     ;
}

