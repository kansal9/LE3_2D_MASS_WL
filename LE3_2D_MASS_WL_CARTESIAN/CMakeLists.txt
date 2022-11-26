CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

#===============================================================================
# Load elements_subdir macro here
# Examples:
#   For declaring a project module:
#         elements_subdir(ElementsExamples)
#===============================================================================
elements_subdir(LE3_2D_MASS_WL_CARTESIAN)

#===============================================================================
# Load elements_depends_on_subdirs macro here
#   For creating a dependency onto an other accessible module
#         elements_depends_on_subdirs(ElementsKernel)
#===============================================================================
elements_depends_on_subdirs(ElementsKernel)
elements_depends_on_subdirs(ST_DM_FilenameProvider)
elements_depends_on_subdirs(LE3_2D_MASS_WL_UTILITIES)
elements_depends_on_subdirs(LE3_2D_MASS_WL_SPHERICAL)

#===============================================================================
# Add the find_package macro (a pure CMake command) here to locate the
# libraries.
# Examples:
#          find_package(CppUnit)
#===============================================================================
find_package(Cfitsio)

#===============================================================================
# Declare the library dependencies here
# Example:
#         elements_add_library(ElementsExamples src/Lib/*.cpp
#                     LINK_LIBRARIES Boost ElementsKernel
#                     INCLUDE_DIRS Boost ElementsKernel
#                     PUBLIC_HEADERS ElementsExamples)
#===============================================================================
elements_add_library(LE3_2D_MASS_WL_CARTESIAN src/lib/*.cpp
                     INCLUDE_DIRS Cfitsio boost ElementsKernel
                     LINK_LIBRARIES ElementsKernel ST_DM_FilenameProvider LE3_2D_MASS_WL_UTILITIES LE3_2D_MASS_WL_SPHERICAL pthread fftw3_threads fftw3 readline
                     PUBLIC_HEADERS LE3_2D_MASS_WL_CARTESIAN)

#===============================================================================
# Declare the executables here
# Example:
# elements_add_executable(ElementsProgramExample src/Program/ProgramExample.cpp
#                        LINK_LIBRARIES Boost ElementsExamples
#                        INCLUDE_DIRS Boost ElementsExamples)
#===============================================================================
elements_add_executable(LE3_2D_MASS_WL_CartesianWL src/program/CartesianWL.cpp
                     INCLUDE_DIRS ElementsKernel LE3_2D_MASS_WL_CARTESIAN
                     LINK_LIBRARIES ElementsKernel LE3_2D_MASS_WL_CARTESIAN)
elements_add_executable(LE3_2D_MASS_WL_PatchesToSphere src/program/PatchesToSphere.cpp
                     INCLUDE_DIRS ElementsKernel LE3_2D_MASS_WL_CARTESIAN
                     LINK_LIBRARIES ElementsKernel LE3_2D_MASS_WL_CARTESIAN)
elements_add_executable(LE3_2D_MASS_WL_Filtering src/program/Filtering.cpp
                     INCLUDE_DIRS ElementsKernel LE3_2D_MASS_WL_CARTESIAN
                     LINK_LIBRARIES ElementsKernel LE3_2D_MASS_WL_CARTESIAN)
#===============================================================================
# Declare the Boost tests here
# Example:
# elements_add_unit_test(BoostClassExample tests/src/Boost/ClassExample_test.cpp
#                       EXECUTABLE BoostClassExample_test
#                       INCLUDE_DIRS ElementsExamples
#                       LINK_LIBRARIES ElementsExamples TYPE Boost)
#===============================================================================
elements_add_unit_test(CartesianParam tests/src/CartesianParam_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_CARTESIAN_CartesianParam_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_CARTESIAN ElementsServices
                     TYPE Boost)
elements_add_unit_test(GenericMap tests/src/GenericMap_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_CARTESIAN_GenericMap_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_CARTESIAN ElementsServices
                     TYPE Boost)
elements_add_unit_test(Projection tests/src/Projection_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_CARTESIAN_Projection_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_CARTESIAN ElementsServices
                     TYPE Boost)
elements_add_unit_test(InpaintingAlgo tests/src/InpaintingAlgo_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_CARTESIAN_InpaintingAlgo_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_CARTESIAN ElementsServices
                     TYPE Boost)
elements_add_unit_test(MapMaker tests/src/MapMaker_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_CARTESIAN_MapMaker_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_CARTESIAN ElementsServices
                     TYPE Boost)
elements_add_unit_test(MatrixProcess tests/src/MatrixProcess_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_CARTESIAN_MatrixProcess_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_CARTESIAN ElementsServices
                     TYPE Boost)
elements_add_unit_test(CartesianAlgoKS tests/src/CartesianAlgoKS_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_CARTESIAN_CartesianAlgoKS_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_CARTESIAN ElementsServices
                     TYPE Boost)
elements_add_unit_test(ShearMap tests/src/ShearMap_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_CARTESIAN_ShearMap_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_CARTESIAN ElementsServices
                     TYPE Boost)
elements_add_unit_test(ConvergenceMap tests/src/ConvergenceMap_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_CARTESIAN_ConvergenceMap_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_CARTESIAN ElementsServices
                     TYPE Boost)
elements_add_unit_test(FilterMR tests/src/FilterMR_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_CARTESIAN_FilterMR_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_CARTESIAN ElementsServices
                     TYPE Boost)
elements_add_unit_test(ReconstructMR tests/src/ReconstructMR_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_CARTESIAN_ReconstructMR_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_CARTESIAN ElementsServices
                     TYPE Boost)
elements_add_unit_test(PatchesToSphereAlgo tests/src/PatchesToSphereAlgo_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_CARTESIAN_PatchesToSphereAlgo_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_CARTESIAN ElementsServices
                     TYPE Boost)
elements_add_unit_test(CartesianProcessor tests/src/CartesianProcessor_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_CARTESIAN_CartesianProcessor_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_CARTESIAN ElementsServices
                     TYPE Boost)
#===============================================================================
# Use the following macro for python modules, scripts and aux files:
#  elements_install_python_modules()
#  elements_install_scripts()
#===============================================================================


#===============================================================================
# Declare the Python programs here
# Examples :
# elements_add_python_program(PythonProgramExample
#                             ElementsExamples.PythonProgramExample)
#===============================================================================

#===============================================================================
# Add the elements_install_conf_files macro
# Examples:
#          elements_install_conf_files()
#===============================================================================
elements_install_conf_files()
