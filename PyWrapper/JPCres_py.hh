#include <memory>
//#include <math.h>
#include <boost/python.hpp>

#include "PwaUtils/DataUtils.hh"

using namespace std;
using namespace boost::python;

class JPCres_py {

public:
  
  JPCres_py();
  ~JPCres_py(){};
  void validLSs();


private:
  std::shared_ptr<jpcRes> _motherRes;
  std::shared_ptr<jpcRes> _daughterRes1;
  std::shared_ptr<jpcRes> _daughterRes2;
};

BOOST_PYTHON_MODULE(JPCres_py)
{
    class_<JPCres_py>("JPCres_py")
      .def("validLSs", &JPCres_py::validLSs)
      ;
}

