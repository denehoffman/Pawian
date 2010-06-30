#include <algorithm>
#include <boost/shared_ptr.hpp>


namespace pawian {
  namespace Collection {

    struct PtrLess {
       template<class PtrType>
      bool operator()(PtrType ptr1, PtrType ptr2) const {
        return (*ptr1) < (*ptr2);
      }
    };

    struct SharedPtrLess {
       template<class PtrType>
      bool operator()(boost::shared_ptr<PtrType> shptr1, boost::shared_ptr<PtrType> shptr2) const {
	PtrType* ptr1=shptr1.get();
	PtrType* ptr2=shptr2.get();
        return (*ptr1) < (*ptr2);
      }
    };

    }

  }
