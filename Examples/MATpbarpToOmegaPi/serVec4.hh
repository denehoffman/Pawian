#ifndef _serVec4_H
#define _serVec4_H

#include <boost/archive/text_oarchive.hpp> 
#include <boost/archive/text_iarchive.hpp>

#include "qft++/tensor/Vector4.hh"

template <typename _Tp> class serVec4 : public Vector4<_Tp> {

  using Vector4<_Tp>::_data;

  friend class boost::serialization::access;

  template<typename Archive>
  void serialize(Archive & ar, const unsigned int) {
    using boost::serialization::make_nvp;

    ar & make_nvp("_data",_data);

    /* Add your own class-variables here in the following way:
      ar & BOOST_SERIALIZATION_NVP(myVar);
      or
      ar & make_nvp("myVar", myVar); // The latter form can be necessary when dealing with templates
    */
  }

  public:
    serVec4() : Vector4<_Tp>::Vector4(){/** Default Constructor */}

  /// Constructor (initialize the 4-vector to be (t,x,y,z))
  serVec4(typename Type<_Tp>::ParamType __t,typename Type<_Tp>::ParamType __x,
	  typename Type<_Tp>::ParamType __y,typename Type<_Tp>::ParamType __z)
    : Vector4<_Tp>::Vector4() {
    this->SetV4(__t,__x,__y,__z);
  }

  /// Copy Constructor
  serVec4(const serVec4<_Tp> &__v4) : Vector4<_Tp>::Vector4(__v4) {}

  /// Copy Constructor
  serVec4(const Vector4<_Tp> &__v4) : Vector4<_Tp>::Vector4(__v4) {}

  virtual ~serVec4(){/** Destructor */}

};

#endif /* _serVec4_H  */

