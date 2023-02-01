#include <vector>

#include "PyWrapper/JPCres_py.hh"
//#include "PwaUtils/DataUtils.hh"

JPCres_py::JPCres_py() :
  _motherRes(std::shared_ptr<jpcRes>(new jpcRes(2,1,1)))
  ,_daughterRes1(std::shared_ptr<jpcRes>(new jpcRes(1,1,-1)))
  ,_daughterRes2(std::shared_ptr<jpcRes>(new jpcRes(1,-1,-1)))
{
}

void JPCres_py::validLSs(){
  std::vector< std::shared_ptr<const LScomb> > theLSVec;
  validLS(_motherRes, _daughterRes1, _daughterRes2, theLSVec);
  std::cout << "for the decay ";
  _motherRes->print(std::cout);
  std::cout  << " -> ";
  _daughterRes1->print(std::cout);
  std::cout << " ";
  _daughterRes2->print(std::cout);
  std::cout << " the valid LS cominations are: " << std::endl;
  std::vector< std::shared_ptr<const LScomb> >::iterator it;
  for (it=theLSVec.begin(); it!=theLSVec.end(); ++it){
    (*it)->print(std::cout);
    std::cout << std::endl;
  }
}


