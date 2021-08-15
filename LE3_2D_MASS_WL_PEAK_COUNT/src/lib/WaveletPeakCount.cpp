/**
 * @file src/lib/WaveletPeakCount.cpp
 * @date 10/07/20
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

#include "LE3_2D_MASS_WL_PEAK_COUNT/WaveletPeakCount.h"
#include <cmath>

using namespace Euclid::WeakLensing::TwoDMass;
using LE3_2D_MASS_WL_CARTESIAN::ShearMap;
using LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap;
static Elements::Logging logger = Elements::Logging::getLogger("WaveletPeakCount");

namespace LE3_2D_MASS_WL_PEAK_COUNT {

WaveletPeakCount::WaveletPeakCount(LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap &convMap,
            LE3_2D_MASS_WL_PEAK_COUNT::PeakParam &PeakParam): m_convMap(convMap),
            m_peakParam(PeakParam), m_MP(convMap.getXdim(), convMap.getYdim()) {
  m_sizeXaxis = convMap.getXdim();
  m_sizeYaxis = convMap.getYdim();
  m_sizeZaxis = convMap.getZdim();
  raMin = convMap.getCoordinateBound().getRaMin();
  raMax = convMap.getCoordinateBound().getRaMax();
  decMin = convMap.getCoordinateBound().getDecMin();
  decMax = convMap.getCoordinateBound().getDecMax();
  m_zMin = convMap.getCoordinateBound().getZMin();
  m_zMax = convMap.getCoordinateBound().getZMax();
}


LE3_2D_MASS_WL_CARTESIAN::Matrix WaveletPeakCount::getSNRimage(LE3_2D_MASS_WL_CARTESIAN::Matrix inputKappaImage,
                                                                 double globalNoise) {
 // divide the signal map by the global noise
 LE3_2D_MASS_WL_CARTESIAN::Matrix mySNRimage = inputKappaImage.multiply(1./globalNoise);
 //logger.info()<<"std dev: "<<globalNoise;
 return mySNRimage;
}

std::vector<std::vector<double> > WaveletPeakCount::getPeaks(std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix> &myBand,
                                                           double globalNoise) {
 //create object to get correct projection
 Projection projection;
 // Get the map characteristics
 double ra0 = 0.5*(raMin + raMax); // in degree
 double dec0 = 0.5*(decMin + decMax); // in degree
 double raRange = (raMax - raMin) * M_PI/180.; //in rad
 double decRange = (decMax - decMin) * M_PI/180.; //in rad

 // Create vectors for peaks data
 std::vector<double> peakX;
 std::vector<double> peakY;
 std::vector<double> peakSNR;
 std::vector<double> min_redshift;
 std::vector<double> max_redshift;
 std::vector<double> scale;
 //logger.info()<<"BandSize: "<<myBand.size();

 logger.info()<<"std dev: "<<globalNoise;
 // Loop on every scale
 for (unsigned int band=0; band<myBand.size()-1; band++) {
  // Divide the image by the noise, to get an SNR image
  LE3_2D_MASS_WL_CARTESIAN::Matrix mySNRimage = getSNRimage(myBand[band],
                        (globalNoise*Euclid::WeakLensing::TwoDMass::norm[band]));
  // Loop over all pixels
  for (unsigned int i=0; i<m_sizeXaxis; i++) {
   for (unsigned int j=0; j<m_sizeYaxis; j++) {
    // Check if any pixel is a local maximum
    if (mySNRimage.isLocalMax(i, j)) {
     //if ((mySNRimage.getValue(i, j)) > m_peakParam.getMinPeakThresh()) {
       // Perform transform from pixel location to ra and dec
       double tmpx = (i+0.5)*raRange/m_convMap.getXdim()-0.5*raRange;
       double tmpy = (j+0.5)*decRange/m_convMap.getYdim()-0.5*decRange;
       std::pair<double, double> radec = projection.getInverseGnomonicProjection(tmpx, tmpy, ra0, dec0);
       //std::pair<double, double> radec = projection.getInverseGnomonicProjection(0.5*raRange*M_PI/180,
                                                                  //0.5*decRange*M_PI/180, ra0, dec0);
       // Save the data into the vectors
  //   logger.info()<<"val: "<<mySNRimage.getValue(0, 0);
     logger.info()<<"PixelSize: "<<m_convMap.getPixelSize();
       double theta = pow(2, band+1)*m_convMap.getPixelSize()*60.;
       peakX.push_back(radec.first);
       peakY.push_back(radec.second);
       peakSNR.push_back(mySNRimage.getValue(i, j));
       scale.push_back(theta);
       min_redshift.push_back(m_zMin);
       max_redshift.push_back(m_zMax);
     //}
    }
   }
  }
 }

 // Create a vector containing all data
 std::vector<std::vector<double> > inputData;

 // And fill this vector before returning it
 inputData.push_back(peakX);//ra
 inputData.push_back(peakY);//dec
 inputData.push_back(peakSNR);
 inputData.push_back(scale);
 inputData.push_back(min_redshift);
 inputData.push_back(max_redshift);

 return inputData;
}

void WaveletPeakCount::savePeakCatalog(const std::string& filename){
  // Allocate memory for kappa and DCT kappa
  LE3_2D_MASS_WL_CARTESIAN::Matrix kappaE(m_sizeXaxis, m_sizeYaxis);
  //LE3_2D_MASS_WL_CARTESIAN::Matrix kappaB(m_sizeXaxis, m_sizeYaxis);

  logger.info()<<"Initializing kappa map";
  //logger.info()<<"binval: " <<m_convMap.getBinValue(1000, 1000, 0);
  for (unsigned int j=0; j<m_sizeYaxis; j++)
  {
    for (unsigned int i=0; i<m_sizeXaxis; i++)
    {
      kappaE.setValue(i, j, m_convMap.getBinValue(i, j, 0));
      //kappaB.setValue(i, j, m_convMap.getBinValue(i, j, 1));
    }
  }
 // Get the noise on the input kappa
 double stdevNoise = (kappaE).getStandardDeviation();
 int nbScales = m_peakParam.getnbScales();
/*
  // Compute the number of scales
  unsigned int nbScales = int(log(m_sizeXaxis)/log(2.))-3.-2.;
  logger.info()<<"number of scales: "<<nbScales=;
*/

  // Perform the kappa map decomposition
  std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix> myBand = m_MP.transformBspline(kappaE, nbScales);

  // Get the peaks data
  std::vector<std::vector<double> > inputData = getPeaks(myBand, stdevNoise);
  logger.info()<<"number of detected peaks: "<<inputData[2].size();
  writePeakCatalog (filename, inputData);
}

void writePeakCatalog (const std::string& filename, std::vector<std::vector<double> > Data) {
  // TODO In next version of EL_FitsIO(3.2) will be available the function to append the columns
 // if (access( filename.c_str(), F_OK ) != -1 ) {
 //   MefFile peakfile(filename, MefFile::Permission::Edit);
 // } else {
    MefFile peakfile(filename, MefFile::Permission::Overwrite);
    const auto col1 = generateColumn<double>("RA_OBJ", Data[0]);
    const auto col2 = generateColumn<double>("DEC_OBJ", Data[1]);
    const auto col3 = generateColumn<double>("ZMIN", Data[4]);
    const auto col4 = generateColumn<double>("ZMAX", Data[5]);
    const auto col5 = generateColumn<double>("THETA", Data[3]);
    const auto col6 = generateColumn<double>("SNR", Data[2]);
  //}
  peakfile.assignBintableExt("PEAK_CATALOG", col1, col2, col3, col4, col5, col6);
}

}  // namespace LE3_2D_MASS_WL_PEAK_COUNT
