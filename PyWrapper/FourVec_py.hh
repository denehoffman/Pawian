#include <iostream>
#include <math.h>
#include <boost/python.hpp>

using namespace std;
using namespace boost::python;

class FourVec_py {
    private:
        double px, py, pz, E;

    public:
        void set (double _px, double _py, double _pz, double _E);
        double mag();
        void greet();
};

BOOST_PYTHON_MODULE(FourVec_py)
{

    class_<FourVec_py>("FourVec_py")
        .def("set", &FourVec_py::set)
        .def("mag", &FourVec_py::mag)
        .def("greet", &FourVec_py::greet)
    ;
}

