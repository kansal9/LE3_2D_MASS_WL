/**
 * @file src/lib/GetMap.cpp
 * @date 05/13/19
 * @author user
 *
 * @copyright (C) 2012-2020 Euclid Science Ground Segment
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 3.0 of the License, or (at your option)
 * any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */

#include "LE3_2D_MASS_WL_CARTESIAN/GetMap.h"

using namespace Euclid::WeakLensing::TwoDMass;
using Euclid::FitsIO::Record;
using Euclid::FitsIO::VecColumn;
using Euclid::FitsIO::VecRaster;
using Euclid::FitsIO::MefFile;
using Euclid::FitsIO::BintableHdu;
using Euclid::FitsIO::ImageHdu;

static Elements::Logging logger = Elements::Logging::getLogger("GetMap");
using namespace Euclid::WeakLensing::TwoDMass;
namespace LE3_2D_MASS_WL_CARTESIAN {

  GetMap::~GetMap() {
   delete m_mapValues;
   m_mapValues = nullptr;
  }

  GetMap::GetMap(double* array, LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam, int nbGalaxies):
                     nGalaxies(nbGalaxies), m_CB(0., 0., 0., 0., 0., 0.), sizeXaxis(1024), m_PixelSize(0.001),
                     sizeYaxis(1024), sizeZaxis(3) {
  sizeXaxis = cartesianParam.getXaxis();
  sizeYaxis = cartesianParam.getYaxis();
  typedef boost::multi_array<double, 3>::index index;

  // Declare the map of data
  m_mapValues = new boost::multi_array<double, 3>(boost::extents[sizeXaxis][sizeYaxis][sizeZaxis]);

  // Assign values from the input array to the map
  for (index k = 0; k != sizeZaxis; ++k){
    for (index j = 0; j != sizeYaxis; ++j){
      for (index i = 0; i != sizeXaxis; ++i){
        (*m_mapValues)[i][j][k] = array[sizeXaxis*sizeYaxis*k+sizeXaxis*j+i];
      }
    }
  }
}

 void GetMap::thresholding(double value){
  typedef boost::multi_array<double, 3>::index index;

  // For each value in the map, add the offset value
    for (index j = 0; j != sizeYaxis; ++j){
      for (index i = 0; i != sizeXaxis; ++i) {
       if (fabs((*m_mapValues)[i][j][0])<value) {
        (*m_mapValues)[i][j][0] = 0;
      }
    }
  }
 }

double GetMap::getSigma() const{
  typedef boost::multi_array<double, 3>::index index;
  double mean(0);
  double squareMean(0);
  unsigned int counts(0);
 // Loop over all values to find the max
  for (index j = 0; j != sizeYaxis; ++j){
   for (index i = 0; i != sizeXaxis; ++i){
    mean += (*m_mapValues)[i][j][0];
    squareMean += (*m_mapValues)[i][j][0]*(*m_mapValues)[i][j][0];
    counts++;
   }
  }
  double stdev = sqrt(squareMean/counts - mean*mean/counts/counts);
  return stdev;
}

 GetMap::GetMap(double* array, int Xaxis, int Yaxis, int Zaxis, int nbGalaxies):
                 sizeXaxis(Xaxis), sizeYaxis(Yaxis), sizeZaxis(Zaxis), nGalaxies(nbGalaxies),
                                  m_PixelSize(0.001), m_CB(0., 0., 0., 0., 0., 0.) {
  typedef boost::multi_array<double, 3>::index index;

  // Declare the map of data
  m_mapValues = new boost::multi_array<double, 3>(boost::extents[sizeXaxis][sizeYaxis][sizeZaxis]);

  // Assign values from the input array to the map
  for (index k = 0; k != sizeZaxis; ++k){
    for (index j = 0; j != sizeYaxis; ++j){
      for (index i = 0; i != sizeXaxis; ++i){
        (*m_mapValues)[i][j][k] = array[sizeXaxis*sizeYaxis*k+sizeXaxis*j+i];
      }
    }
  }
}

 GetMap::GetMap(double* array, int Xaxis, int Yaxis, int Zaxis,
                LE3_2D_MASS_WL_CARTESIAN::CoordinateBound& coordBound, int nbGalaxies):
                 sizeXaxis(Xaxis), sizeYaxis(Yaxis), sizeZaxis(Zaxis), nGalaxies(nbGalaxies),
                 m_CB(coordBound), m_PixelSize(0.001) {
  typedef boost::multi_array<double, 3>::index index;

  // Declare the map of data
  m_mapValues = new boost::multi_array<double, 3>(boost::extents[sizeXaxis][sizeYaxis][sizeZaxis]);

  // Assign values from the input array to the map
  for (index k = 0; k != sizeZaxis; ++k){
    for (index j = 0; j != sizeYaxis; ++j){
      for (index i = 0; i != sizeXaxis; ++i){
        (*m_mapValues)[i][j][k] = array[sizeXaxis*sizeYaxis*k+sizeXaxis*j+i];
      }
    }
  }
 }

 GetMap::GetMap(const std::string& filename):nGalaxies(0), m_CB(0., 0., 0., 0., 0., 0.),
                             m_PixelSize(0.001), sizeXaxis (1024), sizeYaxis(1024), sizeZaxis(3) {
  // Write the values into the object map
   typedef boost::multi_array<double, 3> array_type;
   typedef array_type::index index;

   Euclid::FitsIO::MefFile fitsfile(filename, Euclid::FitsIO::MefFile::Permission::Edit);
   std::vector<std::string> Hdu_names = fitsfile.readHduNames();
   for (size_t i = 0; i < Hdu_names.size(); i++) {
        //logger.info()<<"Name of HDUs in fits file: "<< Hdu_names[i];
        //const auto &ext = fitsfile.accessFirst<Euclid::FitsIO::ImageHdu>("SHEAR_PATCH");
        const auto &ext = fitsfile.accessFirst<Euclid::FitsIO::ImageHdu>(Hdu_names[i]);
        /*const auto records = ext.parseAllRecords<boost::any>();
        sizeXaxis = (records.as<int>("NAXIS1")).value;
        sizeYaxis = (records.as<int>("NAXIS2")).value;
        sizeZaxis = (records.as<int>("NAXIS3")).value;*/
        const auto ramin = ext.parseRecordOr<double>({"RAMIN", 0});
        const auto ramax = ext.parseRecordOr<double>({"RAMAX", 0});
        const auto decmin = ext.parseRecordOr<double>({"DECMIN", 0});
        const auto decmax = ext.parseRecordOr<double>({"DECMAX", 0});
        const auto zmin = ext.parseRecordOr<double>({"ZMIN", 0});
        const auto zmax = ext.parseRecordOr<double>({"ZMAX", 0});
        nGalaxies = ext.parseRecordOr<int>({"NGAL", 0});
        m_PixelSize = ext.parseRecordOr<double>({"PIXSIZE", 0});
        logger.info()<<"Pixel Size found in input Map: "<< m_PixelSize;
        //logger.info()<<"Number of galaxies found in input Map: "<<nGalaxies;
        //logger.info()<<"raminmax: "<<ramin<<"\t"<<ramax << "xaxis: " << sizeXaxis;
        //read Image
        const auto image = ext.readRaster<double, 3>();
        sizeXaxis = image.shape[0];
        sizeYaxis = image.shape[1];
        sizeZaxis = image.shape[2];

        m_CB = CoordinateBound(ramin, ramax, decmin, decmax, zmin, zmax);
        // Declare the map of data
        m_mapValues = new boost::multi_array<double, 3>(boost::extents[sizeXaxis][sizeYaxis][sizeZaxis]);
        // Assign values from the input array to the map
        for (index k = 0; k != sizeZaxis; ++k) { // sizeZaxis = image.length<2>()
          for (index j = 0; j != sizeYaxis; ++j) { //sizeYaxis = image.length<1>()
            for (index i = 0; i != sizeXaxis; ++i){ //sizeXaxis = image.length<0>()
              (*m_mapValues)[i][j][k] = image[{ i, j, k }];
            }
          }
        }
   }
 }

 GetMap::GetMap(GetMap const& copyMap):sizeXaxis(copyMap.sizeXaxis), sizeYaxis(copyMap.sizeYaxis),
      sizeZaxis(copyMap.sizeZaxis), nGalaxies(copyMap.nGalaxies), m_CB(copyMap.m_CB), m_PixelSize(copyMap.m_PixelSize){
  typedef boost::multi_array<double, 3>::index index;
  // Declare the map of data
  m_mapValues = new boost::multi_array<double, 3>(boost::extents[sizeXaxis][sizeYaxis][sizeZaxis]);

  // Assign values from the input array to the map
  for (index k = 0; k != sizeZaxis; ++k){
    for (index j = 0; j != sizeYaxis; ++j){
      for (index i = 0; i != sizeXaxis; ++i) {
        (*m_mapValues)[i][j][k] = (*copyMap.m_mapValues)[i][j][k];
      }
    }
   }
 }

 std::vector<double> GetMap::getMeanValues(){
  // Create the vector of mean values
  std::vector<double> meanValues;
 // Loop over all the values of the map
  typedef boost::multi_array<double, 3>::index index;

  for (index k = 0; k != sizeZaxis; ++k){
    double mean(0.);
    for (index j = 0; j != sizeYaxis; ++j){
      for (index i = 0; i != sizeXaxis; ++i){
        mean += (*m_mapValues)[i][j][k];
      }
    }
    // Scale the sum to the mean
    mean /= sizeXaxis*sizeYaxis;
    // Add that value to the output vector
    meanValues.push_back(mean);
  }
  return meanValues;
 }

 void GetMap::removeOffset(std::vector<double> offset){
  typedef boost::multi_array<double, 3>::index index;

  // For each value in the map, add the offset value
  for (index k = 0; k != sizeZaxis; ++k){
    for (index j = 0; j != sizeYaxis; ++j){
      for (index i = 0; i != sizeXaxis; ++i) {
        (*m_mapValues)[i][j][k] -= offset[k];
      }
    }
  }
 }

 void GetMap::applyGaussianFilter( float sigma){
  applyGaussianFilter(sigma, sigma);
 }

 void GetMap::applyGaussianFilter(float sigmax, float sigmay){
  double fftFactor = 1.0/sizeXaxis/sizeYaxis;
  // Generate a normalized gaussian kernel with sigma
  boost::multi_array<double, 2> gaussianKernel = makeGaussianKernel(sizeXaxis, sizeYaxis, sigmax, sigmay);

  // Create the complex maps
  fftw_complex* kappa_complex = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*sizeXaxis*sizeYaxis);
  fftw_complex* fft_kappa_complex = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*sizeXaxis*sizeYaxis);

  fftw_complex* kernel_complex = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *sizeXaxis*sizeYaxis);
  fftw_complex* fft_kernel_complex = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *sizeXaxis*sizeYaxis);

  fftw_complex* kappaGauss_complex = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*sizeXaxis*sizeYaxis);
  fftw_complex* fft_kappaGauss_complex = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*sizeXaxis*sizeYaxis);

  // Create the plans for transformation
  fftw_plan plan_k_forward = fftw_plan_dft_2d(sizeXaxis, sizeYaxis, kappa_complex,
                                              fft_kappa_complex, FFTW_FORWARD, FFTW_MEASURE);

  fftw_plan plan_kernel_forward = fftw_plan_dft_2d(sizeXaxis, sizeYaxis, kernel_complex,
                                                    fft_kernel_complex, FFTW_FORWARD, FFTW_MEASURE);

  fftw_plan plan_kappaGauss_backward = fftw_plan_dft_2d(sizeXaxis, sizeYaxis, fft_kappaGauss_complex,
                                                        kappaGauss_complex, FFTW_BACKWARD, FFTW_MEASURE);


  // Fill the complex kappa map with convergence map values and kernel values
  for ( int i=0; i<sizeYaxis; i++)
  {
    for ( int j=0; j<sizeXaxis; j++)
    {
      kappa_complex[i*sizeXaxis +j][0] = (*m_mapValues)[i][j][0];
      kappa_complex[i*sizeXaxis +j][1] = 0;//(*m_mapValues)[i][j][1];
      kernel_complex[i*sizeXaxis +j][0] = gaussianKernel[i][j];
      kernel_complex[i*sizeXaxis +j][1] = 0;//gaussianKernel[i][j];
    }
  }

  // Perform the fourier transform of the complex convergence map ans complex kernel
  fftw_execute(plan_k_forward);
  fftw_execute(plan_kernel_forward);


  // Multiply the gaussian kernel and the convergence map in the fourier space
  for ( int i=0; i<sizeXaxis; i++)
  {
    for (int j=0; j<sizeYaxis; j++)
    {
      fft_kappaGauss_complex[j*sizeXaxis +i][0] = fft_kappa_complex[j*sizeXaxis +i][0]*
                                                  fft_kernel_complex[j*sizeXaxis +i][0];
      fft_kappaGauss_complex[j*sizeXaxis +i][1] = fft_kappa_complex[j*sizeXaxis +i][1]*
                                                  fft_kernel_complex[j*sizeXaxis +i][1];
    }
  }

  // Apply backward Fourier transform to get the filtered convergence map
  fftw_execute(plan_kappaGauss_backward);

 // Fill the convergence map
 for (int i=0; i<sizeYaxis; i++) {
  for ( int j=0; j<sizeXaxis; j++) {
   int ii = i<sizeXaxis/2 ? i+sizeXaxis/2 : i-sizeXaxis/2;
   int jj = j<sizeYaxis/2 ? j+sizeYaxis/2 : j-sizeYaxis/2;
   // (*mapValues)[ii][jj][0]=kappaGauss_complex[j*sizeXaxis +i][0]*fftFactor;
   //contents[ii*sizeXaxis +jj] = kappaGauss_complex[i*sizeXaxis +j][0]*fftFactor;
   //mapArray[ ii*sizeXaxis*sizeYaxis + sizeXaxis * jj + 0] = kappaGauss_complex[j*sizeXaxis +i][0]*fftFactor;
   (*m_mapValues)[ii][jj][0]=kappaGauss_complex[j*sizeXaxis +i][0]*fftFactor;
   //logger.info()<<"mapval: "<<(*mapValues)[ii][jj][0]<<"\t"<<"Cont_arr: "<<contents[jj*sizeXaxis +ii];
  //mapArray[ ii*sizeXaxis*sizeYaxis + sizeXaxis * jj + 1] = kappaGauss_complex[j*sizeXaxis +i][1]*fftFactor;
  }
 }

  // free memory
  fftw_destroy_plan(plan_k_forward);
  fftw_destroy_plan(plan_kernel_forward);
  fftw_destroy_plan(plan_kappaGauss_backward);

  fftw_free(kappa_complex);
  fftw_free(fft_kappa_complex);
  fftw_free(kernel_complex);
  fftw_free(fft_kernel_complex);
  fftw_free(kappaGauss_complex);
  fftw_free(fft_kappaGauss_complex);
  fftw_cleanup();
 }

boost::multi_array<double, 2> GetMap::makeGaussianKernel(long sizeX, long sizeY, double sigmaX, double sigmaY) {
 // Create the multi_array output gaussian kernel
 boost::multi_array<double, 2> kernel(boost::extents[sizeX][sizeY]);
 typedef boost::multi_array<double, 2>::index index;
 // Fill in the gaussian kernel
 double sum(0.);
 for (index j = 0; j != sizeY; ++j){
  float y(j-0.5*(sizeY-1));
  for (index i = 0; i != sizeX; ++i) {
   float x(i-0.5*(sizeX-1));
   kernel[i][j]=exp(-(x*x/(2*sigmaX*sigmaX) + y*y/(2*sigmaY*sigmaY)));///(2.*M_PI*sigma*sigma);
   sum+=kernel[i][j];
  }
 }
 double normSum(0.);
 // Normalize the kernel
 for (index j = 0; j != sizeY; ++j) {
  for (index i = 0; i != sizeX; ++i) {
   kernel[i][j] /= sum;
   normSum +=kernel[i][j];
   //logger.info()<<kernel[i][j]<<" ";
  }
 }
//logger.info()<<"normalized sum: "<<normSum;
 return kernel;
}

 int GetMap::getNumberOfGalaxies() const {
  return nGalaxies;
 }

 double GetMap::getPixelSize() const {
  return m_PixelSize;
 }

 double GetMap::getBinValue( int binx, int biny, int binz) const{
   // Make sure the bin values are not wrong
   if (binx >= sizeXaxis) {
    binx = sizeXaxis-1;
   }
   if (biny >= sizeYaxis) {
    biny = sizeYaxis-1;
   }
   if (binz >= sizeZaxis) {
    binz = sizeZaxis-1;
   }
  // Return the value of the map for the (binx, biny, binz)
  return (*m_mapValues)[binx][biny][binz];
}

void GetMap::add_borders(){
 sizeXaxis *= 2;
 sizeYaxis *= 2;
 typedef boost::multi_array<double, 3>::index index;
 // Declare a new array with the right dimensions
 boost::multi_array<double, 3> *borderedMap =
             new boost::multi_array<double, 3>(boost::extents[sizeXaxis][sizeYaxis][sizeZaxis]);
 // Assign values from the old array to the new bordered one
 for (index k = 0; k != sizeZaxis; ++k){
  for (index j = 0; j != sizeYaxis; ++j){
   for (index i = 0; i != sizeXaxis; ++i){
    if (i<sizeXaxis/4 || i>=3*sizeXaxis/4 || j<sizeYaxis/4 || j>=3*sizeYaxis/4){
          (*borderedMap)[i][j][k] = 0.;
    } else {
          (*borderedMap)[i][j][k] = (*m_mapValues)[i-sizeXaxis/4][j-sizeYaxis/4][k];
    }
   }
  }
 }
  // Delete the old map and assign the new one as a member
  delete m_mapValues;
  m_mapValues = borderedMap;
}

void GetMap::remove_borders(){
 sizeXaxis /= 2;
 sizeYaxis /= 2;
 typedef boost::multi_array<double, 3>::index index;
 // Declare a new array with the right dimensions
 boost::multi_array<double, 3> *borderedMap =
             new boost::multi_array<double, 3>(boost::extents[sizeXaxis][sizeYaxis][sizeZaxis]);
 // Assign values from the old array to the new without borders
 for (index k = 0; k != sizeZaxis; ++k) {
  for (index j = 0; j != sizeYaxis; ++j){
   for (index i = 0; i != sizeXaxis; ++i){
     (*borderedMap)[i][j][k] = (*m_mapValues)[i+sizeXaxis/2][j+sizeYaxis/2][k];
   }
  }
 }
 // Delete the old map and assign the new one as a member
 delete m_mapValues;
 m_mapValues = borderedMap;
}

int GetMap::getXdim() const
{
  return sizeXaxis;
}

int GetMap::getYdim() const
{
  return sizeYaxis;
}

int GetMap::getZdim() const
{
  return sizeZaxis;
}

CoordinateBound GetMap::getCoordinateBound() const
{
  return m_CB;
}

unsigned int GetMap::pixelate(float galaxiesPerBin) {
 // Compute the current number of galaxies per bin if known
 float curGalPerBin = float(nGalaxies)/(sizeXaxis*sizeYaxis);
 // Check if the current number of galaxies per bin is lower than the one needed
 // if not, return 1
 if (curGalPerBin>galaxiesPerBin) {
    return 1;
  }
 // Compute here the pixelation....
 // Check what binning would make it closer to the expected value
  bool goodBinning(false);
  int binning(1);
  while (goodBinning==false)
  {
    if (fabs(curGalPerBin*binning*4-galaxiesPerBin)>=fabs(galaxiesPerBin-curGalPerBin))
    {
      goodBinning=true;
    }
    else
    {
      binning *= 4;
    }
  }

  // If the best choice if 1, then do not change anything and return 1
  if (binning == 1)
  {
    return 1;
  }

  // If the best choice means to have one or less pixels then do nothing and return 1
  if (sizeXaxis/int(sqrt(binning)) <= 1 || sizeYaxis/int(sqrt(binning)) <= 1)
  {
    return 1;
  }

 // Create a buffer array to reshape the member array
  boost::multi_array<double, 3> *buffMapValues = new boost::multi_array<double, 3>
    (boost::extents[sizeXaxis/int(sqrt(binning))][sizeYaxis/int(sqrt(binning))][sizeZaxis]);

  typedef boost::multi_array<double, 3>::index index;

  // Reassign binned values to the map
  for (index k = 0; k != sizeZaxis; ++k) {
    for (index j = 0; j != sizeYaxis; ++j) {
      for (index i = 0; i != sizeXaxis; ++i) {
        (*buffMapValues)[i/int(sqrt(binning))][j/int(sqrt(binning))][k] += (*m_mapValues)[i][j][k];
      }
    }
  }

  // Delete the old array
  delete m_mapValues;
  m_mapValues = buffMapValues;

  // Update the axis dimensions of the map
  sizeXaxis /= int(sqrt(binning));
  sizeYaxis /= int(sqrt(binning));

  return int(sqrt(binning));
}

VecRaster<float, 3> GetMap::createRaster() {
  //! [Create and fill a raster]
  VecRaster<float, 3> raster({ sizeXaxis, sizeYaxis, sizeZaxis });
  for (int z = 0; z < sizeZaxis; ++z) {
    for (int y= 0; y < sizeYaxis; ++y) {
      for (int x = 0; x < sizeXaxis; ++x) {
        raster[{ x, y, z }] = (*m_mapValues)[x][y][z];
      }
    }
  }
  return raster;
}

void GetMap::getArray(double *array){
  // Assign values to the elements
  for ( int k = 0; k != sizeZaxis; ++k) {
    for ( int j = 0; j != sizeYaxis; ++j)  {
      for ( int i = 0; i != sizeXaxis; ++i) {
        array[sizeXaxis*sizeYaxis*k+sizeXaxis*j+i] = (*m_mapValues)[i][j][k];
      }
    }
  }
}

bool GetMap::pixelate(int xBinning, int yBinning)
{
  if (xBinning==0 && yBinning==0)
  {
    return false;
  }

  // Put the binning from power of two to actual values
  xBinning = pow(2, xBinning);
  yBinning = pow(2, yBinning);

  // Check the asked binning is not higher or equal to the number of current bins
  if (xBinning>sizeXaxis || yBinning>sizeYaxis)
  {
    return false;
  }

// Create a buffer array to reshape the member array
  boost::multi_array<double, 3> *buffMapValues = new boost::multi_array<double, 3>
    (boost::extents[sizeXaxis/xBinning][sizeYaxis/yBinning][sizeZaxis]);

  typedef boost::multi_array<double, 3>::index index;

  // Reassign binned values to the map
  for (index k = 0; k != sizeZaxis; ++k) {
    for (index j = 0; j != sizeYaxis; ++j){
      for (index i = 0; i != sizeXaxis; ++i) {
        (*buffMapValues)[i/xBinning][j/yBinning][k] += (*m_mapValues)[i][j][k]/xBinning/yBinning;
      }
    }
  }

  // Delete the old array
  delete m_mapValues;
  m_mapValues = buffMapValues;

  sizeXaxis /= xBinning;
  sizeYaxis /= yBinning;

  return true;
}

void GetMap::writeImageHeader(const RecordHdu &hdu, LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam) {

 float cd11 = (m_CB.getRaMax()-m_CB.getRaMin())/sizeXaxis;
 float cd22 = (m_CB.getDecMax()-m_CB.getDecMin())/sizeYaxis;
 float cd12 = 0.0;
 float cd21 = 0.0;
 float crval1 = (floor(sizeXaxis/2)+0.5)*cd11 + m_CB.getRaMin();
 float crval2 = (floor(sizeYaxis/2)+0.5)*cd22 + m_CB.getRaMax();
 float zmin, zmax;
 if (cartesianParam.getParaFileType() == "Conv_Cluster") { //TODO
  zmin = m_CB.getZMin();
  zmax = m_CB.getZMax();
  // image.set("CL_ID", getCL_ID(), "Id of cluster that processed");
  // image.set("CL_MASS", getCL_MASS(), "Mass of the processed cluster");
 }
std::vector<Record<boost::any>> records = {
    { "WCSAXES", 2, "", "Number of axes in World Coordinate System" },
    { "CRPIX1", "", "", "Pixel coordinate of reference point" },
    { "CRPIX2", "", "", "Pixel coordinate of reference point" },
    { "PC1_1", cd11, "", "Coordinate transformation matrix element" },
    { "PC1_2", cd12, "", "Coordinate transformation matrix element" },
    { "PC2_1", cd21, "", "Coordinate transformation matrix element" },
    { "PC2_2", cd22, "", "Coordinate transformation matrix element" },
    { "CDELT1", "", "deg", "Coordinate increment at reference point" },
    { "CDELT2", "", "deg", "Coordinate increment at reference point" },
    { "CUNIT1", "deg", "", "Unit of the first coordinate value" },
    { "CUNIT2", "deg", "", "Unit of the second coordinate value" },
    { "CTYPE1", "RA---TAN", "", "Right ascension, gnomonic projection" },
    { "CTYPE2", "DEC--TAN", "", "Declination, gnomonic projection" },
    { "CRVAL1", crval1, "deg", "Coordinate value at reference point" },
    { "CRVAL2", crval2, "deg", "Coordinate value at reference point" },
    { "LONPOLE", "", "deg", "Native longitude of celestial pole" },
    { "LATPOLE", "", "deg", "Native latitude of celestial pole" },
    { "RADESYS", "", "", "Equatorial coordinate system" },
    { "EQUINOX", "", "", "Equinox of celestial coordinate system (e.g. 2000)" },
    { "DATE-OBS", "", "", "Start of observation" }, //format ‘yyyy-mm-ddThh:mm:ss.sss’
    { "DATE-END", "", "", "End of observation" },
    { "SHE_SOFT", "LENSMC", "", "Origin of Shear Catalog" },
    { "SHE_SEL", "FITCLASS=0", "", "Origin of Shear Catalog" },
    { "SOFTNAME", "LE3_2D_MASS_WL_KS", "", "Software used to create the product" },
    { "SOFTVERS", SWVersion, "", "Software version" },
    { "PROJ", "TAN", "", "projection used" },
    { "PWIDTH", cartesianParam.getPatchWidth(), "", "PatchWidth in degrees" },
    { "PIXSIZE", cartesianParam.getPixelsize(),"", "Pixelsize" },
    { "NITREDSH", cartesianParam.getNItReducedShear(), "", "Number of iterations for reduced shear" },
    { "STDREDSH", cartesianParam.getRSSigmaGauss(),"", "gaussian smoothing sigma for reduced shear" },
    { "NITINP", cartesianParam.getNInpaint(),"", "Number of iterations for inpainting" },
    { "NSCINP", cartesianParam.getnbScales(), "Number of scales for inpainting"},
    { "VARPERSC", cartesianParam.getEqualVarPerScale(), "True if equal variance per scale forced in"},
    { "FBMODE", cartesianParam.getForceBMode(), "" , "True if B-mode forced to zero in the gaps"},
    { "ADDBORD", cartesianParam.get_addBorders(), "" , "True if Borders added"},
    { "DENTYPE", "GAUSSIAN","", "denoising type"},
    { "GAUSSSTD", cartesianParam.getSigmaGauss(), "", "Standard deviation of gaussian smoothing"},
    { "FDRVAL", cartesianParam.getThreshold(),"", "false discovery rate threshold"},
    { "NZBINS", cartesianParam.getnbZBins(),"", "number of redshift bins"},
    { "BAL_BINS", cartesianParam.get_BalancedBins(), "", "True if balanced bins are applied"},
    { "NRESAMPL", cartesianParam.getNSamples(), "", "number of resampling of input dictionary"},
    { "ZMIN", m_CB.getZMin(), "deg", "Minimum Redshift"},
    { "ZMAX", m_CB.getZMax(), "deg", "Maximum Redshift"},
    //{ "NPATCHES", "", "", "number of patches used for the decomposition"},
    //{ "NGAL", int(this->getNumberOfGalaxies()),"", "Number of galaxies"},
  };
  hdu.writeRecords(records);

}

void GetMap::writeMap(const std::string& filename, LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam){
  MefFile f(filename, MefFile::Permission::Overwrite);
  const auto raster = createRaster();
  const auto shape = raster.shape;
  logger.info() << "writing records to primary HDU";
  const auto &primary = f.accessPrimary<>();
  writeImageHeader(primary, cartesianParam);
  logger.info() << "Assigning new Image HDU";
  //const auto &ext = f.initImageExt<float, 3>(cartesianParam.getExtName(), shape);
  const auto &ext = f.assignImageExt(cartesianParam.getExtName(), raster);

 // Image HDU
 //if (cartesianParam.getParaFileType() == "Conv_Patch" || cartesianParam.getExtName() == "SHEAR_PATCH") {
  //Adding extra keys to image header
  ext.writeRecords<double, double>({ "ZMIN", m_CB.getZMin() }, { "ZMAX", m_CB.getZMax() });
 //}
 //if (cartesianParam.getExtName() == "SHEAR_PATCH") {
  ext.writeRecords<double, double>({ "RAMIN", m_CB.getRaMin() }, { "RAMAX", m_CB.getRaMax() });
  ext.writeRecords<double, double> ({ "DECMIN", m_CB.getDecMin() }, { "DECMAX", m_CB.getDecMax() });
  ext.writeRecords<int, double>({ "NGAL",int(this->getNumberOfGalaxies()) },
                                { "PIXSIZE", cartesianParam.getPixelsize() });
// }
}
}  // namespace LE3_2D_MASS_WL_CARTESIAN
