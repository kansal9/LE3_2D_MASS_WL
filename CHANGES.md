# Change log


## 2.6

### New features
* 

### Unknown Bugs
* MRLens Filtering (error between python code and C++ code)

### Expected changes for next release version
* 

### Missing features which will probably be added in next release version
* Use of OU-PHZ inputs (need to coordinate with OU-PHZ people)
* Use mask to estimate number of galaxies per pixel and verify requirements

### Expected validation for next release version
* Validation of Peak Catalogue using full sky Convergence Map
* Validation of reduced shear in Cartesian pipeline
* Validation of MRLens Filtering
* Validation of Monte Carlo maps in Cartesian and Spherical pipeline
* Validation of Shear peak catalogue

### Dependency changes



## 2.5.0


### Work with System Team
* Asked to add few new features to Eulcid library (EL_FitsIO)


### New features
* Extension of the Cartesian convergence to Full sky
* Peak Catalogue built from the full sky convergence map
* MRLens Filtering added
* Input Visibility Mask (from VMPZ_ID) added (read and adapt resolution based on input data)
* New features Patches2Sphere (Merging patches to healpix map)
* Use of SDC-FR to test the code.

### Unknown Bugs
* MRLens Filtering (error between python code and C++ code)

### Maturity Level
* According to current Matrices ~ ML3A

### Missing features which will probably not be added in next release version
* Use of OU-PHZ inputs (need to coordinate with OU-PHZ people)
* Use mask to estimate number of galaxies per pixel and verify requirements

### Expected validation which will probably be done in next release version
* Validation of reduced shear in Cartesian pipeline
* Validation of MRLens Filtering
* Validation of Monte Carlo maps in Cartesian and Spherical pipeline
* Validation of Shear peak catalogue
* Validation of Peak Catalogue using full sky Convergence Map

### Dependency changes
* Update to Elements 5.12.0
* Update to DataModel 8.0.5
* Update to DataModelBindings 8.0.5
* Update to DataModelTools 8.0.5
* Updated to EL_FitsIO 3.1.0

### Added Dependency
* DataSync

### Some tests and modifications needed
* File added at LE3_2D_MASS_WL_UTILITIES/auxdir/*.cpp for PeakCountShear program that needs
  to be changed after next version of EL_FitsIO (Please delete file from auxdir after the changes made)
* Search for TODO inside the project for more things that needs to be done.


## 2.4.0

### After ML2B Review
* Code redesigned to answer ML2B review (not major changes).

### Work with System Team
* Issue opened to add Fits file basic features to Eulcid library (EL_FitsIO)

### New features
* Extension of Cartesian pipeline to the sphere (~50% of the features)
* Optical Cluster catalogue input added in the cluster pipeline
* Dependency to cfitsio/ccfits library deleted
* Worked with system team in order to add 2D-MASS features into Euclid libraries (to be used by anyone)
* Created new project for Pipelines.
* Use of Auxilary Directory for Unit tests


### Maturity Level
* ML2B


### Expected changes for next release version 2.5.0
* Build the Peak Catalogue using full sky convergence map
* Extension of Cartesian pipeline to the sphere (remaining features)
* Merging in Spherical pipeline (per patches)
* Add the Visibility Mask input (from VMPZ_ID)
* Add MRLens Filtering


### Dependency changes
* Update to Elements 5.10.0
* Update to Eden 2.1
* Update to DataModel 8.0.3
* Update to DataModelBindings 8.0.3
* Update to DataModelTools 8.0.3

### Added Dependency
* EL_FitsIO 2.1.0


## 2.2

### Re-Design
* more adapted to the PF requirements.


### New features
* Optimisation of the redshift split of the shear catalogue
* Reduced Shear feature in Cartesian pipeline added
* Produce SNR maps and Monte Carlo realizations
* Extension of the Cartesian pipeline to the sphere (version 1.0 features)
* Optimization of the convergence Peak Catalogue
* Build the Peak Catalogue using Shear catalogue (Mass aperture)
* New project for Verification/Validation created
* Use of Euclid Archieve system to get data and upload output products
* Use of Data Model, Data Model Bindings and Data Model Tools
* use of IAL


### Maturity Level
* ML2B


### Expected changes for next release version 2.4.0
* Build the Peak Catalogue using full sky convergence map
* Extension of the Cartesian pipeline to the sphere (remaining features)
* Merging in Spherical pipeline (per patches)
* Add an input in cluster pipeline (the optical Cluster catalogue) to produce convergence maps 


### Dependency
* Elements 5.2.2
* Eden 2.0


### Added Dependency
* DataModel 2.5.0
* DataModelBindings 2.5.0
* DataModelTools 2.5.0


## 1.0

### Initial features

* Produce patches of E and B-modes convergence maps using Shear Catalogue
* Split the shear catalog into several shear subcatalogs with a given range in ra, dec and z
* Build the Peak Catalog using small patch convergence maps

### Missing features

* Reduced Shear feature in Cartesian pipeline
* SNR maps and Monte Carlo realizations in Cartesian pipeline
* Peak catalogue pipeline from Shear catalogue
* Merge and split in Spherical pipeline (per patches)
* Full-sky convergence maps in Spherical pipeline (full-sky)
* Use of Data Model, Data Model Bindings and Data Model Tools
* Use of IAL

### Dependency
* Elements 5.2.2
* Eden 2.0

### Maturity Level
* ML2A

