#include "PyWrapper/HeliVec_py.hh"
#include "PwaUtils/KinUtils.hh"

HeliVec_py::HeliVec_py() :
  _motherRefQftVec4(Vector4<double>(3.5, 1.0, 0.0, 0.0))
  , _refQftVec4(Vector4<double>(2.8, 0.5, -0.2, -0.1))
  ,_motherQftVec4(Vector4<double>(2.0, 0.3, 0.4, 0.5))
  ,_daughterQftVec4(Vector4<double>(0.5, 0.05, 0.2, 0.1))
{
}

void HeliVec_py::calcHeliVec() {
  std::cout << "_motherRefQftVec4: " << _motherRefQftVec4 << std::endl;
  std::cout << "_refQftVec4: " << _refQftVec4 << std::endl; 
  std::cout << "_motherQftVec4: " << _motherQftVec4 << std::endl;
  std::cout << "_daughterQftVec4: " << _daughterQftVec4 << std::endl;
  Vector4<double> daughterHeliVec4Qft = KinUtils::heliVec(_motherRefQftVec4, _refQftVec4, _motherQftVec4, _daughterQftVec4);
  std::cout << "daughterHeliVec4Qft: " << daughterHeliVec4Qft << std::endl;
}


