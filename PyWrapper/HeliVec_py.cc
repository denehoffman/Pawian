#include "PyWrapper/HeliVec_py.hh"
#include "PwaUtils/KinUtils.hh"

HeliVec_py::HeliVec_py() :
  _motherRefQftVec4(Vector4<double>(3.5, 1.0, 0.0, 0.0))
  , _refQftVec4(Vector4<double>(2.8, 0.5, -0.2, -0.1))
  ,_motherQftVec4(Vector4<double>(2.0, 0.3, 0.4, 0.5))
  ,_daughterQftVec4(Vector4<double>(0.5, 0.05, 0.2, 0.1))
{
}

boost::python::tuple HeliVec_py::calcHeliVec() {
  std::cout << "_motherRefQftVec4: " << _motherRefQftVec4 << std::endl;
  std::cout << "_refQftVec4: " << _refQftVec4 << std::endl; 
  std::cout << "_motherQftVec4: " << _motherQftVec4 << std::endl;
  std::cout << "_daughterQftVec4: " << _daughterQftVec4 << std::endl;
  Vector4<double> daughterHeliVec4Qft = KinUtils::heliVec(_motherRefQftVec4, _refQftVec4, _motherQftVec4, _daughterQftVec4);
  std::cout << "daughterHeliVec4Qft: " << daughterHeliVec4Qft << std::endl;
  return boost::python::make_tuple(daughterHeliVec4Qft.Px(), daughterHeliVec4Qft.Py(), daughterHeliVec4Qft.Pz(), daughterHeliVec4Qft.E());
}


void HeliVec_py::setMother(double px, double py, double pz , double e){
  _motherQftVec4=Vector4<double>(e, px, py, pz);
}

void HeliVec_py::setDaughter(double px, double py, double pz , double e){
  _daughterQftVec4=Vector4<double>(e, px, py, pz);
}


void HeliVec_py::setMotherRef(double px, double py, double pz , double e){
  _motherRefQftVec4=Vector4<double>(e, px, py, pz);
}

void HeliVec_py::setRef(double px, double py, double pz , double e){
  _refQftVec4=Vector4<double>(e, px, py, pz);
}

boost::python::tuple HeliVec_py::getMother(){
  return boost::python::make_tuple(_motherQftVec4.Px(), _motherQftVec4.Py(), _motherQftVec4.Pz(), _motherQftVec4.E());
}

boost::python::tuple HeliVec_py::getDaughter(){
  return boost::python::make_tuple(_daughterQftVec4.Px(), _daughterQftVec4.Py(), _daughterQftVec4.Pz(), _daughterQftVec4.E());
}


boost::python::tuple HeliVec_py::getMotherRef(){
  return boost::python::make_tuple(_motherRefQftVec4.Px(), _motherRefQftVec4.Py(), _motherRefQftVec4.Pz(), _motherRefQftVec4.E());
}

boost::python::tuple HeliVec_py::getRef(){
  return boost::python::make_tuple(_refQftVec4.Px(), _refQftVec4.Py(), _refQftVec4.Pz(), _refQftVec4.E());
}



