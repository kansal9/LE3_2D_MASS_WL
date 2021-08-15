CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

#===============================================================================
# Load elements_subdir macro here
# Examples:
#   For declaring a project module:
#         elements_subdir(ElementsExamples)
#===============================================================================
elements_subdir(LE3_2D_MASS_WL_UTILITIES)

#===============================================================================
# Load elements_depends_on_subdirs macro here
#   For creating a dependency onto an other accessible module
#         elements_depends_on_subdirs(ElementsKernel)
#===============================================================================
elements_depends_on_subdirs(ElementsKernel)
elements_depends_on_subdirs(ST_DataModelBindings)
elements_depends_on_subdirs(ST_DM_HeaderProvider)
elements_depends_on_subdirs(EL_FitsFile)
elements_depends_on_subdirs(EL_FitsData)
elements_depends_on_subdirs(EL_CfitsioWrapper)
#===============================================================================
# Add the find_package macro (a pure CMake command) here to locate the
# libraries.
# Examples:
#          find_package(CppUnit)
find_package(HealpixCpp REQUIRED)
#===============================================================================

#===============================================================================
# Declare the library dependencies here
# Example:
#         elements_add_library(ElementsExamples src/Lib/*.cpp
#                     INCLUDE_DIRS Boost ElementsKernel
#                     LINK_LIBRARIES Boost ElementsKernel
#                     PUBLIC_HEADERS ElementsExamples)
#===============================================================================
elements_add_library(LE3_2D_MASS_WL_UTILITIES src/lib/*.cpp
                     INCLUDE_DIRS ElementsKernel boost HealpixCpp EL_FitsFile ST_DataModelBindings ST_DM_HeaderProvider
                     LINK_LIBRARIES ElementsKernel HealpixCpp EL_FitsFile ST_DataModelBindings ST_DM_HeaderProvider
                     PUBLIC_HEADERS LE3_2D_MASS_WL_UTILITIES)

#===============================================================================
# Declare the executables here
# Example:
# elements_add_executable(ElementsProgramExample src/Program/ProgramExample.cpp
#                        INCLUDE_DIRS Boost ElementsExamples
#                        LINK_LIBRARIES Boost ElementsExamples)
#===============================================================================

#===============================================================================
# Declare the Boost tests here
# Example:
# elements_add_unit_test(BoostClassExample tests/src/Boost/ClassExample_test.cpp
#                       EXECUTABLE BoostClassExample_test
#                       INCLUDE_DIRS ElementsExamples
#                       LINK_LIBRARIES ElementsExamples TYPE Boost)
#===============================================================================
elements_add_unit_test(CatalogData tests/src/CatalogData_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_UTILITIES_CatalogData_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_UTILITIES ElementsServices
                     TYPE Boost)
elements_add_unit_test(DmInput tests/src/DmInput_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_UTILITIES_DmInput_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_UTILITIES ElementsServices
                     TYPE Boost)
elements_add_unit_test(DmOutput tests/src/DmOutput_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_UTILITIES_DmOutput_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_UTILITIES ElementsServices
                     TYPE Boost)
elements_add_unit_test(Utils tests/src/Utils_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_UTILITIES_Utils_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_UTILITIES ElementsServices
                     TYPE Boost)
elements_add_unit_test(NoisyCatalogData tests/src/NoisyCatalogData_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_UTILITIES_NoisyCatalogData_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_UTILITIES
                     TYPE Boost)
elements_add_unit_test(ReadCatalog tests/src/ReadCatalog_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_UTILITIES_ReadCatalog_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_UTILITIES ElementsServices
                     TYPE Boost)
elements_add_unit_test(VisibilityMask tests/src/VisibilityMask_test.cpp 
                     EXECUTABLE LE3_2D_MASS_WL_UTILITIES_VisibilityMask_test
                     LINK_LIBRARIES LE3_2D_MASS_WL_UTILITIES ElementsServices
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
#elements_install_conf_files()
