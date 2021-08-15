LE3-2DMassMapping
 
# Software identification
- Processing function name : 2D-MASS-WL
- Projectname : LE3_2D_MASS_WL_KS
- Profile     : prototype
- Version     : 2.5

# Goal
The processing function (PF) 2D-MASS-WL computes from the shear catalogue:

1. Convergence (E/B mode) maps for small patches
2. Spherical (E/B mode) Convergence maps
3. Peak Catalogue


# Documentation 
| RD | Title | Reference | Issue | Date |
|:--:| :-----| :---------: | :-----: | :----: |
| 1  | ML0 Datapackage    | SEDI-EUCLID-SP-0131-16 | 1.0 | 04/2016 |
| 2  | Mass Mapping output Document  for ML0  | SEDI-EUCLID-SP-0132-16 | 1.0 | 04/2016 |
| 3  | Mass Mapping Requirements | SEDI-EUCLID-SP-0133-16 | 4.0 | 05/2017 |


# Short description
Software 2D-MASS-WL version 2.x which has been redesigned to include new features.
The 2D-MASS-WL PF aims to produce:
-	R-LE3-PRD-F-030: E and B-modes tomographic 2D convergence maps for Nz ≥ 1 redshift bins with a minimum effective angular resolution of 2 arcmin.
    o	Small-field convergence maps 
    o	Full-sky convergence maps
-	G-LE3-PRD-F-031: Related to R-LE3-PRD-F-030, N>100 resamples of the shear catalogue shall be used to compute error bars on the solution
-   R-LE3-PRD-F-032: Catalogues containing the positions and SNRs of peaks detected in 2D Weak Lensing maps (to study peak counts)

See also [RD1] and [RD2] for more details.

# Main programs available
The processing function (PF) 2D-MASS-WL contains four different modules:

+ Cartesian mass mapping module : The processing steps to perform the mass mapping in the plane.
+ Spherical mass mapping module : The processing steps to perform the mass mapping in the sphere.
+ Peak Count module : The processing steps to build the peak catalogue.
+ Utilities: It contains common classes used by all other modules

# How to produce the outputs ?

**E and B modes convergence maps**

In order to produce an E an B modes converegnce maps, we need:
- to split the catalogue into several subcatalogues.
This is done by "CatalogSplitter" module [optional].
- to create shear maps from a shear catalogue.
This is done by the "LE3_2D_MASS_WL_CartesianMapMaker" or "LE3_2D_MASS_WL_SphericalMapMaker". 
- to perform the mass inversion from shear maps to build the convergence maps. 
This is done by the "LE3_2D_MASS_WL_CartesianMassMapping"  or "LE3_2D_MASS_WL_SphericalMassMapping".

**SNR maps**

In order to produce the SNR maps, we need:

- Randomise the shear orientation of each galaxy in the shear catalogue.
- to create noise shear maps from a randomised shear.
This is done by the "LE3_2D_MASS_WL_CartesianMapMaker" or "LE3_2D_MASS_WL_SphericalMapMaker". 
- to perform the mass inversion from noise shear maps to build the noise convergence maps. 
This is done by the "LE3_2D_MASS_WL_CartesianMassMapping"  or "LE3_2D_MASS_WL_SphericalMassMapping".

**Peak catalogue**

Peak catalogue can be produced either using Convergence map or Shear.
In order to produce Peak Catalogue from shear, we need :
- To get Peak catalogue from shear catalog which can be obtained by class "LE3_2D_MASS_WL_PeakCountShear.cpp"


In order to produce Peak Catalogue from convergence, we need :
- To create shear maps from a shear catalogue. This is done by the "LE3_2D_MASS_WL_CartesianMapMaker" or "LE3_2D_MASS_WL_SphericalMapMaker".
- To perform the mass inversion from shear maps to convergence maps. This is done by the "LE3_2D_MASS_WL_CartesianMassMapping" or "LE3_2D_MASS_WL_SphericalMassMapping".
- To get Peak catalogue from convergence map, which can be obtained by either "LE3_2D_MASS_WL_PeakCountConvergence" in case of Cartesian maps or by "LE3_2D_MASS_WL_PeakCountSphere" in case of full sky.


# Configuration

The distribution have been checked with the libraries hereafter:

**'Elements' distribution**

```
- LODEEN 2.1.2
- Elements 5.12.0
```

# Installation with 'Elements'

**Download the project**

```
cd $HOME/Work/Projects
git clone https://gitlab.euclid-sgs.uk/PF-LE3-WL/LE3_2D_MASS_WL_KS.git
cd LE3_2D_MASS_WL_KS
```

**Building the project**

```
make clean
make purge
make configure
make
make install
```

**Running the unitary tests**

```
make test
```

# Data Model

Here is the link to the official data model:
http://euclid.esac.esa.int/dm/dpdd/latest/le3dpd/wl/2D-mass-wl/2D-mass-wlindex.html

**_Input file format_**

We have two different inputs (depend on the pipeline).

First Input: Shear catalog (coming from OU-SHE and OU-PHZ) is expected to contain the following columns:

|**catalog entries**           |                  **Description**                                             | 
| ---------------------------: | :----------------------------------------------------------------------------| 
|RA                            | Right ascension (unit defined as a parameter)                                | 
|DEC                           | Declination (unit defined as a parameter)                                    |
|G1                            | first component of ellipticity (mode, median, mean for instance)             |
|G2                            | second component of ellipticity (mode, median, mean for instance)            | 
|WEIGHT                        | weighting for the galaxy                                                     |
|PHZ_MEDIAN                    | The median of the PHZ_PDF                                                    | 
|PHZ_CORRECTION                | The shift to apply to the PDZs of all objects having this bias identifier    | 
|PHZ_PDF                       | PHZ PDF values, for Z in range [0,6] with 0.01 step                          | 

Note: All columns (except PHZ_MEDIAN & PHZ_PDF) are depend on the different Shear estimation methods (LENSMC, MOMENTSML, KSB or REGAUSS).

Second Input: Cluster catalog (coming from CAT-CL) is expected to contain the following columns:

|**catalog entries**           |                  **Description**                                                | 
| ---------------------------: | :-------------------------------------------------------------------------------| 
|ID                            | Identifier of the cluster candidate                                             | 
|RA                            | Right ascension the cluster center (unit defined as a parameter)                | 
|DEC                           | Declination of the cluster center (unit defined as a parameter)                 |
|Z                             | Estimate of the cluster redshift                                                |
|RADIUS                        | Radius associated to the detection                                              | 
|RICHNESS                      | Richness parameter (a redshift-independent quantity to be used as a mass proxy) |

Third Input: Binary Visibility Mask (in Healpix format coming from VMPZ-ID) is expected to contain the following information:

|**catalog entries**           |                  **Description**                                                | 
| ---------------------------: | :-------------------------------------------------------------------------------| 
|PIXEL                         | Index of the HEALPix pixel cell                                                 | 
|WEIGHT                        | Mask value, ranging from 0 to 1                                                 | 

**_Patch-Based E/B Convergence Maps Parameters_**

|  **Keyword name**  | **Description**                                             | 
| -----------------: | :---------------------------------------------------------- | 
|Nbins               | number of redshift bins                                     | 
|ZMin                | min redshift used                                           | 
|ZMax                | max redshift used                                           |
|BalancedBins        | Repartition of galaxies in redshift bins                    |
|Project             | type of projection                                          | 
|Longitude           | Right ascension of the center of the patch                  | 
|Latitude            | Declination of the center of the patch                      | 
|PixelSize           | Size of the pixels                                          |
|PatchWidth          | Size of the patchs                                          | 
|NItReducedShear     | number of iteration for reduced shear correction            |
|NInpaint            | number of iteration for inpainting (0 if no inpainting)     | 
|NInpScale           | number of wavelet scales for inpainting (0 if no inpainting)| 
|AddBorder           | Adding border to the map (false if not adding border)       | 
|EqualVarPerScale    | Equal Variance inside/outside mask per wavelet scale        | 
|ForceBMode          | flag specifying if B-Mode forced to 0 in the gaps           | 
|GaussSTD            | standard deviation for Gaussian smoothing (0 if no Gaussian smoothing)    | 
|DenoisingAlgo       | denoising type (e.g. Gaussian, name of third party software)| 
|ThresholdFDR        | false discovery rate threshold                              | 


**_Patch-Based SNR Maps Parameters_**

|  **Keyword name**  | **Description**                                             | 
| -----------------: | :---------------------------------------------------------- | 
|Nbins               | number of redshift bins                                     | 
|ZMin                | min redshift used                                           | 
|ZMax                | max redshift used                                           |
|BalancedBins        | Repartition of galaxies in redshift bins                    |
|Project             | type of projection                                          | 
|Longitude           | Right ascension of the center of the patch                  | 
|Latitude            | Declination of the center of the patch                      | 
|PixelSize           | Size of the pixels                                          |
|PatchWidth          | Size of the patchs                                          | 
|NItReducedShear     | number of iteration for reduced shear correction            |
|NInpaint            | number of iteration for inpainting (0 if no inpainting)     | 
|NInpScale           | number of wavelet scales for inpainting (0 if no inpainting)| 
|AddBorder           | Adding border to the map (false if not adding border)       | 
|EqualVarPerScale    | Equal Variance inside/outside mask per wavelet scale        | 
|ForceBMode          | flag specifying if B-Mode forced to 0 in the gaps           | 
|NResamples          | number of resampling of input dictionary                    | 


**_Patch-Based E/B Convergence Maps Parameters for Cluster Pipeline_**

|  **Keyword name**  | **Description**                                             | 
| -----------------: | :---------------------------------------------------------- | 
|MassThreshold       | Select cluster with mass greater than this threshold.       | 
|ZMargin             | Redshift margin for clusters (convergence computed from ZCluster+ZMargin to zmax)  | 
|ZMax                | Maximal redshift considered in the catalog (convergence computed from ZCluster+ZMargin to zmax)  |
|ZMaxHalo            | Maximal redshift of the halo to be kept in the cluster catalog selection  |
|Project             | type of projection                                          | 
|Longitude           | Right ascension of the center of the patch                  | 
|Latitude            | Declination of the center of the patch                      | 
|PixelSize           | Size of the pixels                                          |
|PatchWidth          | Size of the patchs                                          | 
|NItReducedShear     | number of iteration for reduced shear correction            |
|NInpaint            | number of iteration for inpainting (0 if no inpainting)     | 
|NInpScale           | number of wavelet scales for inpainting (0 if no inpainting)| 
|AddBorder           | Adding border to the map (false if not adding border)       | 
|EqualVarPerScale    | Equal Variance inside/outside mask per wavelet scale        | 
|ForceBMode          | flag specifying if B-Mode forced to 0 in the gaps           | 
|GaussSTD            | standard deviation for Gaussian smoothing (0 if no Gaussian smoothing)    | 
|DenoisingAlgo       | denoising type (e.g. Gaussian, name of third party software)| 
|ThresholdFDR        | false discovery rate threshold                              | 


**_Patch-Based SNR Maps Parameters for Cluster Pipeline_**

|  **Keyword name**  | **Description**                                             | 
| -----------------: | :---------------------------------------------------------- | 
|MassThreshold       | Select cluster with mass greater than this threshold.       | 
|ZMargin             | Redshift margin for clusters (convergence computed from ZCluster+ZMargin to zmax)   | 
|ZMax                | Maximal redshift considered in the catalog (convergence computed from ZCluster+ZMargin to zmax)   |
|ZMaxHalo            | Maximal redshift of the halo to be kept in the cluster catalog selection  |
|Project             | type of projection                                          | 
|Longitude           | Right ascension of the center of the patch                  | 
|Latitude            | Declination of the center of the patch                      | 
|PixelSize           | Size of the pixels                                          |
|PatchWidth          | Size of the patchs                                          | 
|NItReducedShear     | number of iteration for reduced shear correction            |
|NInpaint            | number of iteration for inpainting (0 if no inpainting)     | 
|NInpScale           | number of wavelet scales for inpainting (0 if no inpainting)| 
|AddBorder           | Adding border to the map (false if not adding border)       | 
|EqualVarPerScale    | Equal Variance inside/outside mask per wavelet scale        | 
|ForceBMode          | flag specifying if B-Mode forced to 0 in the gaps           | 
|NResamples          | number of resampling of input dictionary                    | 


**_Spherical E/B Convergence Maps Parameters_**

|  **Keyword name**  | **Description**                                             | 
| -----------------: | :---------------------------------------------------------- | 
|Nbins               | number of redshift bins                                     | 
|ZMin                | min redshift used                                           | 
|ZMax                | max redshift used                                           |
|BalancedBins        | Repartition of galaxies in redshift bins                    |
|Nside               | Healpix pixel size                                          |
|PatchWidth          | Size of the patchs                                          | 
|NItReducedShear     | number of iteration for reduced shear correction            |
|NInpaint            | number of iteration for inpainting (0 if no inpainting)     | 
|NInpScale           | number of wavelet scales for inpainting (0 if no inpainting)| 
|AddBorder           | Adding border to the map (false if not adding border)       | 
|EqualVarPerScale    | Equal Variance inside/outside mask per wavelet scale        | 
|ForceBMode          | flag specifying if B-Mode forced to 0 in the gaps           | 
|GaussSTD            | standard deviation for Gaussian smoothing (0 if no Gaussian smoothing)    | 
|DenoisingAlgo       | denoising type (e.g. Gaussian, name of third party software)| 


**_Spherical SNR Maps Parameters_**

|  **Keyword name**  | **Description**                                             | 
| -----------------: | :---------------------------------------------------------- | 
|Nbins               | number of redshift bins                                     | 
|ZMin                | min redshift used                                           | 
|ZMax                | max redshift used                                           |
|BalancedBins        | Repartition of galaxies in redshift bins                    |
|Nside               | Healpix pixel size                                          |
|PatchWidth          | Size of the patchs                                          | 
|NItReducedShear     | number of iteration for reduced shear correction            |
|NInpaint            | number of iteration for inpainting (0 if no inpainting)     | 
|NInpScale           | number of wavelet scales for inpainting (0 if no inpainting)| 
|AddBorder           | Adding border to the map (false if not adding border)       | 
|EqualVarPerScale    | Equal Variance inside/outside mask per wavelet scale        | 
|ForceBMode          | flag specifying if B-Mode forced to 0 in the gaps           | 
|NResamples          | number of resampling of input dictionary                    | 


**_Peak Count using Wavelet Filter Configuration file_**

|  **Keyword name**  | **Description**                                             | 
| -----------------: | :---------------------------------------------------------- | 
|Nbins               | number of redshift bins                                     | 
|ZMin                | min redshift used                                           | 
|ZMax                | max redshift used                                           |
|BalancedBins        | Repartition of galaxies in redshift bins                    |
|Project             | type of projection                                          | 
|Longitude           | Right ascension of the center of the patch                  | 
|Latitude            | Declination of the center of the patch                      | 
|PixelSize           | Size of the pixels                                          |
|PatchWidth          | Size of the patchs                                          | 
|NItReducedShear     | number of iteration for reduced shear correction            |
|NInpaint            | number of iteration for inpainting (0 if no inpainting)     | 
|NInpScale           | number of wavelet scales for inpainting (0 if no inpainting)| 
|AddBorder           | Adding border to the map (false if not adding border)       | 
|EqualVarPerScale    | Equal Variance inside/outside mask per wavelet scale        | 
|ForceBMode          | flag specifying if B-Mode forced to 0 in the gaps           | 
|GaussSTD            | standard deviation for Gaussian smoothing (0 if no Gaussian smoothing)    | 
|DenoisingAlgo       | denoising type (e.g. Gaussian, name of third party software)| 
|ThresholdFDR        | false discovery rate threshold                              | 
|NScale              | Number of scales for the peaks                              | 


**_Peak Count using Mass Aperture Configuration file_**

|  **Keyword name**  | **Description**                                             | 
| -----------------: | :---------------------------------------------------------- | 
|Nbins               | number of redshift bins                                     | 
|ZMin                | min redshift used                                           | 
|ZMax                | max redshift used                                           |
|BalancedBins        | Repartition of galaxies in redshift bins                    |
|Project             | type of projection                                          | 
|Longitude           | Right ascension of the center of the patch                  | 
|Latitude            | Declination of the center of the patch                      | 
|PixelSize           | Size of the pixels                                          |
|PatchWidth          | Size of the patchs                                          | 
|NPeakScale          | List of aperture radius for peaks                           | 


**_Convergence maps Output file_**

- Patches: N 3 dimensional cartesian maps (per redshift bin) containing 2D E and B mode convergence 
maps of the N patches of the sky.

The image contains an image with the number of galaxies per pixel in the maps (can serve to derive a noise level or to identify zones where inpainting has been performed- zone with no galaxy measured).


- Sphere: For each redshift bin (1 hdu per redshift bin), the E/B modes are stored in a 2d binary table with two vector columns (KAPPA_E, KAPPA_B)s.

For each redshift bin, a GALCOUNT column will contain the number of Galaxies measured per pixel on the sphere and redshift bin considered, which can also be used to identify the area inpainted (pixel without measured galaxies).


**_SNR maps Output file_**

- Patches: N 3 dimensional cartesian maps (per redshift bin) containing the 2D SNR associated 
to the E and B modes convergence maps for the N patches of the sky.

- Sphere: For each redshift bin (1 hdu per redshift bin), the maps are stored in a 2d binary table with two vector columns (SNR_E, SNR_B), associated to the E/B SNR maps.


**_Peak Count Output file_**

The output peak catalog contains the following columns:

| **Catalogue entries** |           **Description**     |  
|---------------------:| :--------------------------- | 
|RA_OBJ                 | right ascension of the peak   | 
|DEC_OBJ                | declination of the peak       |
|ZMIN                   | minimum redshift of the bin where the peak has been detected| 
|ZMAX                   | maimum redshift of the bin where the peak has been detected|
|THETA                  | scale where peak has been detected|
|SNR                    | SNR associated to the peak        | 


# Quality checking

This step compiles the code with quality flags and runs many profiling tools. 
A dashboard allow us to have a look on the analysis.

The command to start the analysis locally in LODEEN is: 
checkQuality -m all

The dashboard produced by 'checkQuality' in CODEEN is visible with the url: 

https://codeen.euclid-ec.org/jenkins/job/LE3_2D_MASS_WL_KS/job/develop/

The Maturity level and sonar report is visible with the url:

https://codeen-app.euclid-ec.org/sonar/dashboard?id=LE3_2D_MASS_WL_KS_2.5


# Software available

## C++ binary components
The code should be launched with the standard Euclid procedure with the command:

E-Run <Project> [version] <Executable> options

The LE3_2D_MASS_WL_KS has 4 modules:

o	The Cartesian mass mapping module includes the processing steps to perform the mass mapping in the plane.

o	The Spherical mass mapping module includes the processing steps to perform the mass mapping in the sphere.

o	The Peak Count module includes the processing steps to build the peak catalogue

o	And a Utilities module on top of it.

The executables are based on the operations.

Here is the list of the binary softwares available in the built :

### Cartesian Mass Mapping Module

**I) LE3_2D_MASS_WL_CatalogSplitter**

The module CatalogSplitter extract sub-catalogues from the shear catalogue according to the requirements:

- The redshift range (redshift min and max to be considered)
- The number of bins (in the defined redshift range)
- A Boolean to set the galaxy repartition in the redshift bins (regular redshift bins or with an equal number of bins)

```
NAME
	LE3_2D_MASS_WL_CatalogSplitter --  is used to cut out the catalogues in several subcatalogues
    
SYNOPSIS 

E-Run LE3_2D_MASS_WL_KS 2.5 LE3_2D_MASS_WL_CatalogSplitter

The options available (listed via the --help) are:
	  --workdir=<path>          name of the working directory 
	  --inputCatalog=<name>          input Catalog in fits/xml format
  	  --paramFile=<path_and_name_of_parameter_file>          Input Parameter File to split input catalog
  	  --Sub_Catalogs=<name>          jason/txt file name to store Sub-Catalog Names

```
**Running the code**

**_Split the catalogue into subcatalogues :_**
```
Example:
E-Run LE3_2D_MASS_WL_KS 2.5 LE3_2D_MASS_WL_CatalogSplitter --workdir=/run/media/user/Backup_Drive/WL_Results
--paramFile=ParamConvergencePatch.xml --inputCatalog=euc-test-le3-wl-InputLE2Catalog.xml
```

**II) LE3_2D_MASS_WL_CartesianMapMaker**

The module LE3_2D_MASS_WL_CartesianMapMaker:

- build the shear maps from the shear catalogue. It includes gnomonic projection, projection effects correction and pixelisation.

```
NAME
	LE3_2D_MASS_WL_CartesianMapMaker
    
SYNOPSIS
E-Run LE3_2D_MASS_WL_KS 2.5 LE3_2D_MASS_WL_CartesianMapMaker

The options available (listed via the --help) are:
	  --workdir=<path> 					                	name of the working directory 
	  --input_ShearCatalog=<name>				            input shear Catalog in fits/xml format
      --input_ClusterCatalog=<name>                         input cluster Catalog in fits/xml
	  --paramFile=<path_and_name_of_parameter_file>       	Input Parameter File in XML
	  --outShearMap=<name>   				                output shear Map in fits format
      --outShearMapList=<name>                              List of output shear Maps in text/json format

```

**Running the code**

**_Build the shear maps from the shear catalogue:_**
```
Example:
E-Run LE3_2D_MASS_WL_KS 2.5 LE3_2D_MASS_WL_CartesianMapMaker --workdir=/run/media/user/Backup_Drive/WL_Results
--paramFile=ParamConvergencePatch.xml --input_ShearCatalog=euc-test-le3-wl-InputLE2Catalog.xml
```

**_Build the shear maps from the shear catalogue & cluster catalogue:_**
```
Example:
E-Run LE3_2D_MASS_WL_KS 2.5 LE3_2D_MASS_WL_CartesianMapMaker --workdir=/run/media/user/Backup_Drive/WL_Results
--paramFile=ParamConvergenceCluster.xml --input_ShearCatalog=euc-test-le3-wl-InputLE2Catalog.xml --input_ClusterCatalog=euc-test-le3-wl-twodmass-ClusterCatalog.xml
```

**III) LE3_2D_MASS_WL_CartesianMassMapping**

The LE3_2D_MASS_WL_CartesianMassMapping performs the (Cartesian) direct mass inversion using KS or KS+ methods. The optional input parameters (for KS+):
-	Number of inpainting iterations 
-	Number of reduced shear iterations
-	A Boolean to set B-mode constraint 
-	A Boolean to force the power spectrum constraint to enforce the variance inside the gaps to be equal to the variance outside the gaps at different scales 
-	Number of scales used by the power spectrum constraint 
-	The width of the Gaussian filter for reduced shear iteration
-	A Boolean to add borders to the map to deal with border effects 


```

NAME
	LE3_2D_MASS_WL_CartesianMassMapping
    
SYNOPSIS 
E-Run LE3_2D_MASS_WL_KS 2.5 LE3_2D_MASS_WL_CartesianMassMapping

The options available (listed via the --help) are:
	  --workdir=<path> 						                name of the working directory 
	  --inputShearMap=<name>				            	input Shear Map in fits/json file format
      --inputConvMap=<name>                                 input convergence Map in fits format
	  --paramFile=<path_and_name_of_parameter_file>       	Input Parameter File in XML
      --outConvMap=<name>                                   output convergence Map in fits format
	  --outConvMapXML=<name>   		     	                output convergence Map in xml format
      --outShearMap=<name>                                  output shear Map in fits format

```
**Running the code**

**_Direct Kaiser and Squires, Inverse Kaiser and Squires or Direct Kaiser and Squires including inpainting (based on input parameters in Parameter file) :_**
```
Example:

E-Run LE3_2D_MASS_WL_KS 2.5 LE3_2D_MASS_WL_CartesianMassMapping --workdir=/run/media/user/Backup_Drive/WL_Results
--inputShearMap=ShearMap.fits --paramFile=ParamConvergencePatch.xml

```


### Spherical Mass Mapping Module

**I) LE3_2D_MASS_WL_CatalogSplitter**

The module CatalogSplitter extract a sub-catalogue from the shear catalogue according to the requirements:

- The redshift range (redshift min and max to be considered)
- The number of bins (in the defined redshift range)
- A Boolean to set the galaxy repartition in the redshift bins (regular redshift bins or with an equal number of bins)

```
NAME
	LE3_2D_MASS_WL_CatalogSplitter --  is used to cut out the catalogues in several subcatalogs
    
SYNOPSIS 

E-Run LE3_2D_MASS_WL_KS 2.5 LE3_2D_MASS_WL_CatalogSplitter

The options available (listed via the --help) are:
	  --workdir=<path> 						              name of the working directory 
	  --inputCatalog=<name>					              input Catalog in fits/xml format
  	  --paramFile=<path_and_name_of_parameter_file>       Input Parameter File to split input catalog
  	  --Sub_Catalogs=<name> 			                  jason/txt file name to store Sub-Catalog Names

```
**Running the code**

**_Split the catalogue into subcatalogues :_**
```
Example:
E-Run LE3_2D_MASS_WL_KS 2.5 LE3_2D_MASS_WL_CatalogSplitter --workdir=/run/media/user/Backup_Drive/WL_Results
--paramFile=ParamsConvergenceSphere.xml --inputCatalog=euc-test-le3-wl-InputLE2Catalog.xml
```

**II) LE3_2D_MASS_WL_SphericalMapMaker**

The module LE3_2D_MASS_WL_SphericalMapMaker:

- build the shear maps from the shear catalogue. It includes projection and pixelisation.

```
NAME
	LE3_2D_MASS_WL_SphericalMapMaker
    
SYNOPSIS
E-Run LE3_2D_MASS_WL_KS 2.5 LE3_2D_MASS_WL_SphericalMapMaker

The options available (listed via the --help) are:
	  --workdir=<path> 					                	            name of the working directory 
	  --inputShearCatalog=<name>				                        input Catalog in fits/xml format
	  --sphericalParameterFile=<path_and_name_of_parameter_file>       	Input Parameter File in XML
	  --outShearMap=<name>   				  	                        output shear Map E & B mode in fits format
	  --GalCountMap=<name>   				   	                        output Galaxy count Map which conatins number of Galaxies per pixel for eachredshift bin in XML format

```

**Running the code**

**_Build the shear maps from the shear catalogue:_**
```
Example:
E-Run LE3_2D_MASS_WL_KS 2.5 LE3_2D_MASS_WL_SphericalMapMaker --workdir=/run/media/user/Backup_Drive/WL_Results
--sphericalParameterFile=ParamConvergenceSphere.xml --inputShearCatalog=euc-test-le3-wl-InputLE2Catalog.xml
```

**III) LE3_2D_MASS_WL_SphericalMassMapping**

The LE3_2D_MASS_WL_SphericalMassMapping performs the direct mass inversion using KS or KS+ methods(not implemented yet). The optional input parameters (for KS+):
-	Number of inpainting iterations 
-	Number of reduced shear iterations
-	A Boolean to set B-mode constraint 
-	A Boolean to force the power spectrum constraint to enforce the variance inside the gaps to be equal to the variance outside the gaps at different scales 
-	Number of scales used by the power spectrum constraint 
-	The width of the Gaussian filter for reduced shear iteration
-	A Boolean to add borders to the map to deal with border effects 

```

NAME
	LE3_2D_MASS_WL_SphericalMassMapping
    
SYNOPSIS 
E-Run LE3_2D_MASS_WL_KS 2.5 LE3_2D_MASS_WL_SphericalMassMapping

The options available (listed via the --help) are:
	  --workdir=<path> 				                      	name of the working directory 
	  --inGamma=<name>			                    	    input E/B mode Shear Map in fits(BinTable) format
	  --inKappa=<name>				                        input E/B mode Convergence Map in fits(BinTable) format
	  --paramFile=<path_and_name_of_parameter_file>         Input Parameter File in XML
	  --outputKappa=<name>   				                Output convergence map E & B mode in BinTable Format
	  --outputShear=<name>   			            	    Output shear map E & B mode in BinTable Format
      --outputKappaXML=<name>                               Output XML product for Spherical convergence map

```
**Running the code**

**_Direct Kaiser and Squires, Inverse Kaiser and Squires or Direct Kaiser and Squires including inpainting (based on input parameters in Parameter file) :_**
```
Example:

E-Run LE3_2D_MASS_WL_KS 2.5 LE3_2D_MASS_WL_SphericalMassMapping --workdir=/run/media/user/Backup_Drive/WL_Results
--inGamma=ShearMap.fits --paramFile=ParamConvergenceSphere.xml

```

### Spherical Mass Mapping Module

**I) PeakCount**

This module is supposed to build the peak catalogue.

A peak corresponds to a maximum with the four neighbours lower.
This module compute peaks using two methods:
- from Convergence map and using wavelet filters.
- from the shear catalogue and using mass aperture filters.

```
PeakCount from Convergence Map and using Wavelet Filters

NAME
	LE3_2D_MASS_WL_PeakCountConvergence -- Build the peak catalogue from a convergence map
    
SYNOPSIS 
E-Run LE3_2D_MASS_WL_KS 2.5 LE3_2D_MASS_WL_PeakCountConvergence

 The options available (listed via the --help) are:
   	  --workdir=<path> 				                            	name of the working directory 
	  --inputConvMap=<name>				                            input Convergence Map	
	  --paramPeakConvergence=<path_and_name_of_parameter_file>      input Parameter File in XML
      --outputPConvCatalog=<name>                                   Output Peak Count Catalog in fits
	  --outputPConvCatalogXML=<name>		                    	output Peak Catalog in XML format

```
**Running the code**

**_Peak Catalog estimated from a convergence map :_**
```
Example:
E-Run LE3_2D_MASS_WL_KS 2.5 LE3_2D_MASS_WL_PeakCountConvergence --workdir=/run/media/user/Backup_Drive/WL_Results
--inputConvMap=convergenceMap.fits --paramPeakConvergence=ParamsPeakCountConvergence.xml

```

```
PeakCount from Full sky Convergence Map and using Wavelet Filters

NAME
	LE3_2D_MASS_WL_PeakCountSphere -- Build the peak catalogue from full sky convergence map
    
SYNOPSIS 
E-Run LE3_2D_MASS_WL_KS 2.5 LE3_2D_MASS_WL_PeakCountSphere

 The options available (listed via the --help) are:
   	  --workdir=<path> 				                            	name of the working directory 
	  --inputConvMap=<name>				                            input Convergence Map	
	  --paramPeakConvergence=<path_and_name_of_parameter_file>      input Parameter File in XML
      --outputPeakCatalog=<name>                                   Output Peak Count Catalog in fits
	  --outputPeakCatalogXML=<name>		                    	  output Peak Catalog in XML format

```
**Running the code**

**_Peak Catalog estimated from full sky convergence map :_**
```
Example:
E-Run LE3_2D_MASS_WL_KS 2.5 LE3_2D_MASS_WL_PeakCountSphere --workdir=/run/media/user/Backup_Drive/WL_Results
--inputConvMap=FullSkyConvergenceMap.fits --paramPeakConvergence=ParamsPeakCountConvergence.xml

```

```
PeakCount from Shear Catalog and using Mass Aperture Filters

NAME
	LE3_2D_MASS_WL_PeakCountShear -- Build the peak catalogue from a Shear Catalog
    
SYNOPSIS 
E-Run LE3_2D_MASS_WL_KS 2.5 LE3_2D_MASS_WL_PeakCountShear

 The options available (listed via the --help) are:
 	 --workdir=<path> 				                            	name of the working directory 
	 --inputShearCatalog=<name>			                          	input shear Catalog in fits/xml format
	 --paramPeakMassAperture=<path_and_name_of_parameter_file>      input Parameter File in XML
	 --paramPatchMap=<path_and_name_of_parameter_file>              input Patch Parameter File in XML
     --outputPeakMACatalog=<name>                                   Output Peak Count Catalog in fits
	 --outputPeakMACatalogXML=<name>			                    output Peak Catalog in XML format

```
**Running the code**

**_Peak Catalog estimated from a Shear catalog :_**
```
Example:
E-Run LE3_2D_MASS_WL_KS 2.5 LE3_2D_MASS_WL_PeakCountShear --workdir=/run/media/user/Backup_Drive/WL_Results
--inputShearCatalog=euc-test-le3-wl-InputLE2Catalog.xml --paramPatchMap=ParamConvergencePatch.xml
--paramPeakMassAperture=ParamsPeakCountMassAperture.xml

```
