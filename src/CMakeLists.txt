set(SINGLECYRSTAL_SOURCES posResTest.cxx HagridCrystal.cxx SingleCrystalConstruction.cxx PosResTestActionInitialization.cxx PencilBeamRunAction.cxx PencilBeamRun.cxx PencilBeamEventAction.cxx PencilBeamPrimaryGenerator.cxx PhysicsList.cxx OpticalDetection.cxx)
#set(SINGLECYRSTAL_SOURCES posResTest.cxx FlatSurfaceCurrent.cxx OpticalDetection.cxx SingleCrystalConstruction.cxx PosResTestActionInitialization.cxx PencilBeamRunAction.cxx PencilBeamRun.cxx PencilBeamEventAction.cxx PencilBeamSteppingAction.cxx PencilBeamPrimaryGenerator.cxx PhysicsList.cxx)

add_executable(posResTest ${SINGLECYRSTAL_SOURCES})
target_link_libraries(posResTest ${Geant4_LIBRARIES})
install(TARGETS posResTest DESTINATION ${CMAKE_INSTALL_PREFIX})
