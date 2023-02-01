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
  boost::python::tuple calcHeliVec();
  void setMother(double px, double py, double pz , double e);
  void setDaughter(double px, double py, double pz , double e);
  void setMotherRef(double px, double py, double pz , double e);
  void setRef(double px, double py, double pz , double e);
  boost::python::tuple getMother();
  boost::python::tuple getDaughter();
  boost::python::tuple getMotherRef();
  boost::python::tuple getRef();

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
      .def("setMother", &HeliVec_py::setMother)
      .def("setDaughter", &HeliVec_py::setDaughter)
      .def("setMotherRef", &HeliVec_py::setMotherRef)
      .def("setRef", &HeliVec_py::setRef)
      .def("getMother", &HeliVec_py::getMother)
      .def("getDaughter", &HeliVec_py::getDaughter)
      .def("getMotherRef", &HeliVec_py::getMotherRef)
      .def("getRef", &HeliVec_py::getRef)
      ;
}

