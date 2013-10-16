#pragma once

#include <algorithm>
#include <memory>

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
      bool operator()(const std::shared_ptr<const PtrType>& shptr1, 
		      const std::shared_ptr<const PtrType>& shptr2) const  {
	const PtrType* ptr1=shptr1.get();
	const PtrType* ptr2=shptr2.get();
        return (*ptr1) < (*ptr2);
      }

      template<class PtrType>
      bool operator()(std::shared_ptr<PtrType> shptr1, 
		      std::shared_ptr<PtrType> shptr2) const {
	PtrType* ptr1=shptr1.get();
	PtrType* ptr2=shptr2.get();
        return (*ptr1) < (*ptr2);
      }
    };

  }
}
