CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

#===============================================================================
# Load macros and functions for Elements-based projects
#===============================================================================

find_package(ElementsProject)

#set(ELEMENTS_PARALLEL "" ON)
#remove_definitions(-D_GLIBCXX_PARALLEL)

set(ELEMENTS_CPP14 OFF)
set(ELEMENTS_CPP17 ON
    CACHE STRING "Enable C++17 compilation"
    FORCE)

#===============================================================================
# Declare project name and version
# Example with dependency:
#                         elements_project(MyProject 1.0 USE Element 3.9)
#===============================================================================

elements_project(LE3_2D_MASS_WL_KS 2.6 USE ST_DataModelTools 9.0.3 EleFits 5.0)