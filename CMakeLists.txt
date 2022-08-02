CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

#===============================================================================
# Load macros and functions for Elements-based projects
#===============================================================================

find_package(ElementsProject)

option(ELEMENTS_PARALLEL "" ON)
remove_definitions(-D_GLIBCXX_PARALLEL)

#===============================================================================
# Declare project name and version
# Example with dependency:
#                         elements_project(MyProject 1.0 USE Element 3.9)
#===============================================================================

elements_project(LE3_2D_MASS_WL_KS 2.6 USE ST_DataModelTools 9.0.2 EleFits 5.0)