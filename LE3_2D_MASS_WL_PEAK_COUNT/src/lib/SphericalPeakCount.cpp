/**
 * @file src/lib/SphericalPeakCount.cpp
 * @date 06/18/21
 * @author vanshika kansal
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

#include "LE3_2D_MASS_WL_PEAK_COUNT/SphericalPeakCount.h"

static Elements::Logging logger = Elements::Logging::getLogger("Spherical PeakCount");

using namespace Euclid::WeakLensing::TwoDMass;
using namespace Euclid::WeakLensing::TwoDMass::Spherical;

namespace LE3_2D_MASS_WL_PEAK_COUNT {

SphericalPeakCount::SphericalPeakCount(Healpix_Map<double>& convergenceE, 
                              LE3_2D_MASS_WL_PEAK_COUNT::PeakParam &PeakParam):
      m_peakParam(PeakParam), KappaE(convergenceE) {
   m_nside = KappaE.Nside();
   m_nlmax = 3*m_nside-1;
   m_npix = KappaE.Npix();
   m_order = KappaE.Order();
   logger.info()<<"lmax: "<<m_nlmax;
   logger.info()<<"npix: "<<m_npix;
   m_nbScales = m_peakParam.getnbScales();
}

Healpix_Map<double> SphericalPeakCount::getSNRimage_hp(Healpix_Map<double>& inputKappa,
                                                                 double globalNoise) {
   Healpix_Map<double> mySNRimage;
   mySNRimage.SetNside(m_nside, RING);
   mySNRimage.fill(0.);
 // divide the signal map by the global noise
   for (int i =0; i<m_npix; i++) {
     mySNRimage[i] = inputKappa[i] * (1./globalNoise);
   } 
 return mySNRimage;
}

std::vector<std::vector<double> > SphericalPeakCount::getPeaks_hp(std::vector<Healpix_Map<double> >& myBand,
                                                           double globalNoise) {
 std::vector<double> peakX;
 std::vector<double> peakY;
 std::vector<double> peakSNR;
 std::vector<double> min_redshift;
 std::vector<double> max_redshift;
 std::vector<double> scale;

 logger.info()<<"std dev: "<<globalNoise;

 for (unsigned int band=0; band<myBand.size()-1; band++) {
  // Divide the image by the noise, to get an SNR image
  Healpix_Map<double> mySNRimage = getSNRimage_hp(myBand[band],
                        (globalNoise*Euclid::WeakLensing::TwoDMass::SphericalNorm[band]));
   fix_arr< int, 8 > neighbours;
   double maxVal;
   for (int ipix=0; ipix<m_npix; ipix++) {
     mySNRimage.neighbors(ipix, neighbours);
/*     logger.info()<<"neighbours size: "<< neighbours.size();
     logger.info()<<"neighbours 1: "<< neighbours[0];
     logger.info()<<"neighbours 2: "<< neighbours[1];
     logger.info()<<"neighbours 3: "<< neighbours[2];
     logger.info()<<"neighbours 4: "<< neighbours[3];
     logger.info()<<"neighbours 5: "<< neighbours[4];
     logger.info()<<"neighbours 6: "<< neighbours[5];
     logger.info()<<"neighbours 7: "<< neighbours[6];
     logger.info()<<"neighbours 8: "<< neighbours[7];*/
     std::pair<double, double> radec;
     for (size_t i = 0; i < neighbours.size(); i++) {
        if(mySNRimage[ipix] < neighbours[i]) {
           maxVal = mySNRimage[neighbours[i]];
           radec = getIndex2RaDec(mySNRimage, ipix);
        }
     }

     peakX.push_back(radec.first);
     peakY.push_back(radec.second);
     peakSNR.push_back(maxVal);
     scale.push_back(band);
     min_redshift.push_back(0.);
     max_redshift.push_back(0.);
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

void SphericalPeakCount::savePeakCatalog_hp(const std::string& filename){

  double stdev = getStd(KappaE);

  std::vector<Healpix_Map<double> > band = transformBspline_hp(KappaE, m_nbScales);

  // Get the peaks data
  std::vector<std::vector<double> > inputData = getPeaks_hp(band, stdev);
  logger.info()<<"number of detected peaks: "<<inputData[2].size();
  writePeakCatalog (filename, inputData);
}

}  // namespace LE3_2D_MASS_WL_PEAK_COUNT
