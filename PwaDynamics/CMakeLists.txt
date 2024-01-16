add_library(PwaDynamics  
  BBUnstableParFit.cc  
  BBUnstableParFitFcn.cc
  BarrierFactor.cc
  BreitWignerFunction.cc
  FVector.cc
  FVectorOmnes.cc
  Flatte.cc
  Johnson.cc
  KMatrixBase.cc
  KMatrixFunctions.cc
  KMatrixKPiSFocus.cc
  KMatrixNonRel.cc
  KMatrixOmnes.cc
  KMatrixPiPiS.cc
  KMatrixPiPiSsimple4piPhp.cc
  KMatrixRel.cc
  KMatrixRelBg.cc
  KMatrixSlowAdlerCorRel.cc
  KPole.cc
  KPoleBarrier.cc
  KPoleOmnes.cc
  LASS.cc  
  PPole.cc
  PPoleBarrier.cc
  PPoleOmnes.cc
  PVectorBgOmnes.cc
  PVectorKPiSFocus.cc
  PVectorRel.cc
  PVectorRelBg.cc
  PVectorSlowCorRel.cc
  PhaseSpace4Pi.cc
  PhaseSpaceFactory.cc
  PhaseSpaceIsobar.cc
  PhaseSpaceIsobarAS.cc
  PhaseSpaceIsobarDudek.cc
  PhaseSpaceIsobarDudekAngularMomentum.cc
  PhaseSpaceIsobarDudekUnstableRhoPi.cc
  PhaseSpaceIsobarLUT.cc
  PhaseSpaceIsobarLUT1D.cc
  PhaseSpaceIsobarReid.cc
  PhaseSpaceIsobarReidAngularMomentum.cc
  PhaseSpaceList.cc
  RadMultipoleFormFactor.cc
  TMatrixBase.cc
  TMatrixNonRel.cc
  TMatrixRel.cc
  Voigtian.cc
)

target_include_directories(PwaDynamics
  PRIVATE 
    ${CMAKE_SOURCE_DIR}    
)

target_link_libraries(PwaDynamics  
  Utils
  Particle
  qft++
  qft++Extension
  ErrLogger
  ROOT::Minuit2
)

install(TARGETS PwaDynamics
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

add_executable(BBUnstableParFitApp
  BBUnstableParFitApp.cc
)

target_include_directories(BBUnstableParFitApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}    
)

target_link_libraries(BBUnstableParFitApp
  PRIVATE
    PwaDynamics
    Utils
    Particle
    qft++
    qft++Extension
    ErrLogger
)

add_executable(LUTFileApp
  LUTFileApp.cc
)

target_include_directories(LUTFileApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}    
)

target_link_libraries(LUTFileApp
  PRIVATE
    PwaDynamics
    Utils
    Particle
    qft++
    qft++Extension
    ErrLogger
)

install(TARGETS BBUnstableParFitApp LUTFileApp
  RUNTIME DESTINATION bin  
)