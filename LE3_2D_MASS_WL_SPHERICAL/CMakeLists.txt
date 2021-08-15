CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

#===============================================================================
# Load elements_subdir macro here
# Examples:
#   For declaring a project module:
#         elements_subdir(ElementsExamples)
#===============================================================================
elements_subdir(LE3_2D_MASS_WL_SPHERICAL)

#===============================================================================
# Load elements_depends_on_subdirs macro here
#   For creating a dependency onto an other accessible module
#         elements_depends_on_subdirs(ElementsKernel)
#===============================================================================
elements_depends_on_subdirs(ElementsKernel)
elements_depends_on_subdirs(LE3_2D_MASS_WL_UTILITIES)
#===============================================================================
# Add the find_package macro (a pure CMake command) here to locate the
# libraries.
# Examples:
#          find_package(CppUnit)
find_package(Cfitsio)
find_package(HealpixCpp REQUIRED)
#===============================================================================

#===============================================================================
# Declare the library dependencies here
# Example:
#         elements_add_library(ElementsExamples src/Lib/*.cpp
#                     LINK_LIBRARIES Boost ElementsKernel
#                     INCLUDE_DIRS Boost ElementsKernel
#                     PUBLIC_HEADERS ElementsExamples)
#===============================================================================
elements_add_library(LE3_2D_MASS_WL_SPHERICAL src/lib/*.cpp
                     INCLUDE_DIRS boost HealpixCpp Cfitsio ElementsKernel
                     LINK_LIBRARIES ElementsKernel LE3_2D_MASS_WL_UTILITIES HealpixCpp Cfitsio
                     PUBLIC_HEADERS LE3_2D_MASS_WL_SPHERICAL)

#===============================================================================
# Declare the executables here
# Example:
# elements_add_executable(ElementsProgramExample src/Program/ProgramExample.cpp
#                        LINK_LIBRARIES Boost ElementsExamples
#                        INCLUDE_DIRS Boost ElementsExamples)
#===============================================================================
elements_add_executable(LE3_2D_MASS_WL_SphericalMapMaker src/program/LE3_2D_MASS_WL_SphericalMapMaker.cpp
                     INCLUDE_DIRS ElementsKernel LE3_2D_MASS_WL_SPHERICAL
                     LINK_LIBRARIES ElementsKernel LE3_2D_MASS_WL_SPHERICAL)
elements_add_executable(LE3_2D_MASS_WL_SphericalMassMapping src/program/LE3_2D_MASS_WL_SphericalMassMapping.cpp
                     INCLUDE_DIRS ElementsKernel LE3_2D_MASS_WL_SPHERICAL
                     LINK_LIBRARIES ElementsKernel LE3_2D_MASS_WL_SPHERICAL)
elements_add_executable(LE3_2D_MASS_WL_SphericalMCMaps src/program/LE3_2D_MASS_WL_SphericalMCMaps.cpp
                     INCLUDE_DIRS ElementsKernel LE3_2D_MASS_WL_SPHERICAL
                     LINK_LIBRARIES ElementsKernel LE3_2D_MASS_WL_SPHERICAL)

#===============================================================================
# Declare the Boost tests here
# Example:
# elements_add_unit_test(BoostClassExample tests/src/Boost/ClassExample_test.cpp
#                       EXECUTABLE BoostClassExample_test
#                       INCLUDE_DIRS ElementsExamples
#                       LINK_LIBRARIES ElementsExamples TYPE Boost)
#===============================================================================
elements_add_unit_test(Sph_map_maker tests/src/Sph_map_maker_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_SPHERICAL_Sph_map_maker_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_SPHERICAL
                     TYPE Boost)
elements_add_unit_test(Sph_mass_mapping tests/src/Sph_mass_mapping_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_SPHERICAL_Sph_mass_mapping_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_SPHERICAL
                     TYPE Boost)
elements_add_unit_test(SphericalParam tests/src/SphericalParam_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_SPHERICAL_SphericalParam_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_SPHERICAL
                     TYPE Boost)
elements_add_unit_test(GetSphericalMCMaps tests/src/GetSphericalMCMaps_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_SPHERICAL_GetSphericalMCMaps_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_SPHERICAL
                     TYPE Boost)
elements_add_unit_test(SphericalIO tests/src/SphericalIO_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_SPHERICAL_SphericalIO_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_SPHERICAL
                     TYPE Boost)
elements_add_unit_test(SphericalInpainting tests/src/SphericalInpainting_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_SPHERICAL_SphericalInpainting_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_SPHERICAL ElementsServices
                     TYPE Boost)
elements_add_unit_test(SphericalUtils tests/src/SphericalUtils_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_SPHERICAL_SphericalUtils_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_SPHERICAL ElementsServices
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
