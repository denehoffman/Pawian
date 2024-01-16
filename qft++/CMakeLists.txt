add_library(qft++
  relativistic-quantum-mechanics/BlattWeisskopf.cc
  relativistic-quantum-mechanics/DiracAntiSpinor.cc
  relativistic-quantum-mechanics/DiracSpinor.cc
  relativistic-quantum-mechanics/PolVector.cc
  relativistic-quantum-mechanics/Utils.cc
  tensor/TensorIndex.cc
  tensor/TensorIndexOrder.cc  
)

target_include_directories(qft++
  PRIVATE 
    ${CMAKE_SOURCE_DIR}/  
)

install(TARGETS qft++
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
