/**
 * @file src/lib/MassAperturePeakCount.cpp
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

#include "LE3_2D_MASS_WL_PEAK_COUNT/MassAperturePeakCount.h"
#include <cmath>

using namespace Euclid::WeakLensing::TwoDMass;
using LE3_2D_MASS_WL_CARTESIAN::ShearMap;
using LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap;
static Elements::Logging logger = Elements::Logging::getLogger("MassAperturePeakCount");

namespace LE3_2D_MASS_WL_PEAK_COUNT {

/*MassAperturePeakCount::MassAperturePeakCount(LE3_2D_MASS_WL_CARTESIAN::ShearMap &shearMap,
      LE3_2D_MASS_WL_CARTESIAN::ShearMap &noisedShearMap, LE3_2D_MASS_WL_PEAK_COUNT::PeakParam &PeakParam,
      LE3_2D_MASS_WL_CARTESIAN::CartesianParam& CParam): m_shearMap(shearMap), m_noisedShearMap(noisedShearMap),
       m_peakParam(PeakParam), m_cartesianParam(CParam)*/
MassAperturePeakCount::MassAperturePeakCount(LE3_2D_MASS_WL_CARTESIAN::ShearMap &shearMap,
      LE3_2D_MASS_WL_CARTESIAN::ShearMap &noisedShearMap, LE3_2D_MASS_WL_PEAK_COUNT::PeakParam &PeakParam):
      m_shearMap(shearMap), m_noisedShearMap(noisedShearMap), m_peakParam(PeakParam)
{
  m_sizeXaxis = shearMap.getXdim();
  m_sizeYaxis = shearMap.getYdim();
  m_sizeZaxis = shearMap.getZdim();
  m_raMin = shearMap.getCoordinateBound().getRaMin();
  m_raMax = shearMap.getCoordinateBound().getRaMax();
  m_decMin = shearMap.getCoordinateBound().getDecMin();
  m_decMax = shearMap.getCoordinateBound().getDecMax();
  m_zMin = shearMap.getCoordinateBound().getZMin();
  m_zMax = shearMap.getCoordinateBound().getZMax();
  radius =  m_peakParam.getApPeakRadius();
}

std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix> MassAperturePeakCount::getMassApertureMap
                                  (LE3_2D_MASS_WL_CARTESIAN::ShearMap &inShearMap) {
//GetMap MassAperturePeakCount::getMassApertureMap(LE3_2D_MASS_WL_CARTESIAN::ShearMap &shearMap, double radius){
  LE3_2D_MASS_WL_CARTESIAN::ShearMap m_localShearMap(inShearMap);

  // Allocate memory for shear
  LE3_2D_MASS_WL_CARTESIAN::Matrix shearE(m_sizeXaxis, m_sizeYaxis);
  LE3_2D_MASS_WL_CARTESIAN::Matrix shearB(m_sizeXaxis, m_sizeYaxis);

  //logger.info()<<"binval: " <<m_shearMap.getBinValue(1000, 1000, 0);
  for (unsigned int j=0; j<m_sizeYaxis; j++)
  {
    for (unsigned int i=0; i<m_sizeXaxis; i++)
    {
      shearE.setValue(i, j, m_localShearMap.getBinValue(i, j, 0));
      shearB.setValue(i, j, m_localShearMap.getBinValue(i, j, 1));
    }
  }
  std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix> myBand;
  for (size_t iter = 0; iter<radius.size(); iter++) {
   LE3_2D_MASS_WL_CARTESIAN::Matrix Map(m_sizeXaxis, m_sizeYaxis);
   Map = createMassApertureMap (shearE, shearB, iter);

  /* if (radius[iter] == 4) {
     double *ShearArray = new double[m_sizeXaxis*m_sizeYaxis*m_sizeZaxis];
     for (size_t i=0; i<m_sizeXaxis; i++){
      for (size_t j=0; j<m_sizeYaxis; j++) {
       ShearArray[j*m_sizeXaxis +i] = shearE.getValue(i, j);
       ShearArray[m_sizeXaxis*m_sizeYaxis + j*m_sizeXaxis +i] = shearB.getValue(i, j);
      }
     }

     ShearMap SMap(ShearArray, m_sizeXaxis, m_sizeYaxis, m_sizeZaxis);
    const std::string shMap = "MassApertureMap_Scale_" + std::to_string(radius[iter]) + ".fits";
  // Free the memory
  delete [] ShearArray;
  ShearArray = nullptr;
  SMap.writeMap(shMap, m_cartesianParam);
   }// end if*/
   myBand.push_back(Map);
  }

 return myBand;
}

LE3_2D_MASS_WL_CARTESIAN::Matrix MassAperturePeakCount::createMassApertureMap(LE3_2D_MASS_WL_CARTESIAN::Matrix &shearE,
                                        LE3_2D_MASS_WL_CARTESIAN::Matrix &shearB, unsigned int iter) {
   double rmax = radius[iter]*4.;

   LE3_2D_MASS_WL_CARTESIAN::Matrix Map(m_sizeXaxis, m_sizeYaxis);

   for (int i0 = rmax; i0 < m_sizeXaxis-rmax; i0++) {
    for (int j0 = rmax; j0 < m_sizeYaxis-rmax; j0++) {
      for (int i = 0; i < 2.*rmax; i++) {
        for (int j = 0; j < 2.*rmax; j++) {
           int ii, jj, a, b;
           double a1, a2, phi, kernel, r;
           ii = i0 + (i-rmax);
           jj = j0 + (j-rmax);
           a = (ii-i0);
           b = (jj-j0);
           phi = atan2(a, b);
           r = sqrt(pow(a, 2.)+ pow(b, 2.));
           kernel = 0;
           if (r < rmax) {
             kernel = pow(r, 2)/(4.* M_PI *pow(radius[iter], 4))*exp(-pow(r, 2.)/(2.*pow(radius[iter], 2.)));
           }
           a1 = shearE.getValue(ii, jj)*cos(2.*phi);
           a2 = shearB.getValue(ii, jj)*sin(2.*phi);
           Map.setValue(i0, j0, Map.getValue(i0, j0) + (a1 - a2)*kernel);
           }
         }
      }
    }
   return Map;
}

void MassAperturePeakCount::saveMAPeakCatalog(const std::string& filename) {
//create object to get correct projection
 Projection projection;
 // Get the map characteristics
 double ra0 = 0.5*(m_raMin + m_raMax); // in degree
 double dec0 = 0.5*(m_decMin + m_decMax); // in degree
 double raRange = (m_raMax - m_raMin) * M_PI/180.; //in rad
 double decRange = (m_decMax - m_decMin) * M_PI/180.; //in rad

 // Create vectors for peaks data
 std::vector<double> peakRA;
 std::vector<double> peakDEC;
 std::vector<double> peakSNR;
 std::vector<double> min_redshift;
 std::vector<double> max_redshift;
 std::vector<double> scale;

  // get Mass aperture maps wrt each radius
  std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix> myShearBand = getMassApertureMap(m_shearMap);

  // get Mass aperture maps wrt each radius
  std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix> myNoisedShearBand = getMassApertureMap(m_noisedShearMap);

  for (size_t iter = 0; iter!=myNoisedShearBand.size(); ++iter) {
     double stdevNoise = (myNoisedShearBand[iter]).getStandardDeviation();
     LE3_2D_MASS_WL_CARTESIAN::Matrix mySNRimage = myShearBand[iter].multiply(1./stdevNoise);
     // Loop over all pixels
     for (unsigned int i=0; i<m_sizeXaxis; i++) {
       for (unsigned int j=0; j<m_sizeYaxis; j++) {
        // Check if any pixel is a local maximum
        if (mySNRimage.isLocalMax(i, j)) {
         //if ((mySNRimage.getValue(i, j)) > m_peakParam.getMinPeakThresh()) {
          // Perform transform from pixel location to ra and dec
          double tmpx = (i+0.5)*raRange/m_sizeXaxis-0.5*raRange;
          double tmpy = (j+0.5)*decRange/m_sizeYaxis-0.5*decRange;
          std::pair<double, double> radec = projection.getInverseGnomonicProjection(tmpx, tmpy, ra0, dec0);
          peakRA.push_back(radec.first);
          peakDEC.push_back(radec.second);
          peakSNR.push_back(mySNRimage.getValue(i, j));
//          scale.push_back(1 + int (log(radius[iter])/log(2))); //TODO save values in theta (arcmin)
          scale.push_back(iter);
          min_redshift.push_back(m_zMin);
          max_redshift.push_back(m_zMax);
         //}
        }
       }
     }
  }
  logger.info()<<"number of detected peaks: "<<peakRA.size();
  std::vector<std::vector<double> > Data;

  Data.push_back(peakRA);
  Data.push_back(peakDEC);
  Data.push_back(peakSNR);
  Data.push_back(scale);
  Data.push_back(min_redshift);
  Data.push_back(max_redshift);

  writePeakCatalog (filename, Data);
}

}  // namespace LE3_2D_MASS_WL_PEAK_COUNT
