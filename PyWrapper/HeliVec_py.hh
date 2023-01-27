#include <iostream>
#include <math.h>
#include <boost/python.hpp>

using namespace std;
using namespace boost::python;

struct Vec4_py{
  double px, py, pz, E;

  void printOut(){
    std::cout << "(E,px,py,pz) = ("
	      << E << ","
	      << px << ","
	      << py << ","      
	      << pz << ")"
	      << std::endl;
  }

};

class HeliVec_py {

public:
  
  HeliVec_py();
  ~HeliVec_py(){};
  void calcHeliVec();

private:
  Vec4_py motherRefVec4;
  Vec4_py refVec4;
  Vec4_py motherVec4;
  Vec4_py daughterVec4;

};

BOOST_PYTHON_MODULE(HeliVec_py)
{
    class_<HeliVec_py>("HeliVec_py")
        .def("calcHeliVec", &HeliVec_py::calcHeliVec)
     ;
}

