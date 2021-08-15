CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

#===============================================================================
# Load macros and functions for Elements-based projects
#===============================================================================

find_package(ElementsProject)

option(ELEMENTS_PARALLEL "" ON)

#===============================================================================
# Declare project name and version
# Example with dependency:
#                         elements_project(MyProject 1.0 USE Element 3.9)
#===============================================================================

#elements_project(LE3_2D_MASS_WL_KS 2.3 USE Elements 5.2.2 ST_DataModelBindings 2.5.0 ST_DataModelTools 2.5.0)
#elements_project(LE3_2D_MASS_WL_KS 2.5 USE ST_DataModelTools 8.0.3 EL_FitsIO 2.1.0)
#elements_project(LE3_2D_MASS_WL_KS 2.5 USE ST_DataModelTools 8.0.4 EL_FitsIO 3.0)
#elements_project(LE3_2D_MASS_WL_KS 2.5 USE ST_DataModelTools 8.0.4 EL_FitsIO 2.2.0)
#elements_project(LE3_2D_MASS_WL_KS 2.5 USE ST_DataModelTools 8.0.5 EL_FitsIO 2.2.0)
elements_project(LE3_2D_MASS_WL_KS 2.5 USE ST_DataModelTools 8.0.5 EL_FitsIO 3.1.0)
