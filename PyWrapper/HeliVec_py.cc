// classes example
#include "PyWrapper/HeliVec_py.hh"
#include "qft++/tensor/Vector4.hh"
#include "PwaUtils/KinUtils.hh"

HeliVec_py::HeliVec_py(){
  motherVec4.px=0.3;
  motherVec4.py=0.4;
  motherVec4.pz=0.5;
  motherVec4.E=2.;
  daughterVec4.px=0.05;
  daughterVec4.py=0.2;
  daughterVec4.pz=0.1;
  daughterVec4.E=0.5;
  refVec4.px=0.5;
  refVec4.py=-0.2;
  refVec4.pz=-0.1;
  refVec4.E=2.8;
  motherRefVec4.px=0.0;
  motherRefVec4.py=0.0;
  motherRefVec4.pz=1.0;
  motherRefVec4.E=3.5;
}

void HeliVec_py::calcHeliVec() {
  std::cout << "motherRefVec4: ";
  motherRefVec4.printOut();
  std::cout << "\nrefVec4: ";
  refVec4.printOut();
  std::cout << "\nmotherVec4: ";
  motherVec4.printOut();
  std::cout << "\ndaughterVec4: ";
  daughterVec4.printOut();
  Vector4<double> motherRefQft(motherRefVec4.E, motherRefVec4.px, motherRefVec4.py, motherRefVec4.pz);
  Vector4<double> refVec4Qft(refVec4.E, refVec4.px, refVec4.py, refVec4.pz);
  Vector4<double> daughterVec4Qft(daughterVec4.E, daughterVec4.px, daughterVec4.py, daughterVec4.pz);
  Vector4<double> motherVec4Qft(motherVec4.E, motherVec4.px, motherVec4.py, motherVec4.pz);
  Vector4<double> daughterHeliVec4Qft = KinUtils::heliVec(motherRefQft, refVec4Qft, motherVec4Qft, daughterVec4Qft);
  std::cout << "daughterHeliVec4Qft: " << daughterHeliVec4Qft << std::endl;
}


