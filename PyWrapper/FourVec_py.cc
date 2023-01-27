// classes example
#include "PyWrapper/FourVec_py.hh"
#include "qft++/tensor/Vector4.hh"

void FourVec_py::set (double _px, double _py, double _pz, double _E) {
    px = _px;
    py = _py;
    pz = _pz;
    E = _E;
    return;
}

double FourVec_py::mag() {
  Vector4<double> the4Vec(E,px,py,pz);
  std::cout << "the4Vec mag: " << the4Vec.P() << std::endl;  
  return sqrt(px*px+py*py+pz*pz);
}

void FourVec_py::greet() {
    std::cout << "hello" << std::endl;
}

