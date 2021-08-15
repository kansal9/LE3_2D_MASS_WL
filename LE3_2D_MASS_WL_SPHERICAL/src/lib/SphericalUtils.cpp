/**
 * @file src/lib/SphericalUtils.cpp
 * @date 06/17/21
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

#include "LE3_2D_MASS_WL_SPHERICAL/SphericalUtils.h"

static Elements::Logging logger = Elements::Logging::getLogger("SphericalUtils");

namespace Euclid {
 namespace WeakLensing {
  namespace TwoDMass {
   namespace Spherical {

void readSphericalParameterFile (const boost::filesystem::path& ParamFile,
                                   LE3_2D_MASS_WL_SPHERICAL::SphericalParam &params) {
   logger.info() << "Using Paramfile: " << ParamFile << " as DM input product";
   //check file is of XML format
   if (true == checkFileType(ParamFile.native(), Euclid::WeakLensing::TwoDMass::signXML)) {
    if (true == fileHasField(ParamFile.native(), "DpdTwoDMassParamsConvergenceSphere")) {
      logger.info()<<"Parameter file is for Spherical Convergence..";
      params.getConvergenceSphereParam(ParamFile.native());
      logger.info()<< "Done reading Spherical Convergence Parameter file. . . ";
    }
   } else {
       logger.info()<<"Parameter file is not provided or not in XML format . .";
   }
}

void applyGaussianFilter_hp(Healpix_Map<double>& map, double sigma) {
  int nside = map.Nside();
  arr<double> weight;
  weight.alloc(3*nside-1);
  weight.fill(1.);

  int nlmax = 3*nside-1;

  Alm<xcomplex<double> > map_lm(nlmax, nlmax);
  map_lm.SetToZero();

 // TODO Change sigma from pixel to radian
  double sigma2fwhm = sigma * SIGMA2FWHM;
 // double fwhm2sigma= fwhm / SIGMA2FWHM;

  logger.info()<<"sigma2fwhm: " << sigma2fwhm;
  //logger.info()<<"fwhm2sigma: " << fwhm2sigma;

  map2alm_iter (map, map_lm, 3, weight);

  smoothWithGauss( map_lm, sigma2fwhm); //fwhm needs to be in radian

  map.fill(0.);

  alm2map(map_lm, map);
}

 double b3_spline (double val) {
   double x1 = fabs(pow((val-2), 3));
   double x2 = fabs (pow((val-1), 3));
   double x3 = fabs(pow(val, 3));
   double x4 = fabs(pow((val + 1), 3));
   double x5 = fabs (pow((val + 2), 3));
   double OutValue = 1./12. * (x1 - 4. * x2 + 6. * x3 - 4. * x4 + x5);
   return OutValue;
 }

 void applyThreshold (Alm<xcomplex<double> >& alphaT, double threshVal, int m_nlmax) {
   for (int l =0; l<=m_nlmax; ++l) {
    for (int m = 0; m<=l; ++m){
      if (fabs(alphaT(l, m)) < threshVal) {
         alphaT(l, m) = 0.;
        // alphaT(l, m).real(0.);
        // alphaT(l, m).imag(0.);
      }
    }
   }
 }

 double max_abs_alm(Alm<xcomplex<double> >& alm, int m_nlmax) {
    double max=0.;
    for (int l=0; l <= m_nlmax; ++l) {	
      for (int m=0; m <= l; ++m) {
         if (fabs (real(alm (l,m))) > max) {
            max = fabs (real(alm (l,m)));
         }
         if (fabs (imag(alm (l,m))) > max) {
            max = fabs (imag(alm (l,m)));
         }
      }
    }
	return max;
 }

 double getStd(Healpix_Map<double>& map) {
  double mean(0);
  double squareMean(0);
  unsigned int counts(0);
  int m_npix = map.Npix();
 // Loop over all values to find the max
  for (int i=1; i<m_npix; i++) {
    mean += map[i];
    squareMean += map[i]*map[i];
    counts++;
  }
  double stdev = sqrt(squareMean/counts - mean*mean/counts/counts);

  return stdev;
 }

 void getFilter(double* array, double lc, int m_nlmax) {
   for (int l=0; l<=m_nlmax; l++) {
          array[l]=(3./2.) * b3_spline(2.0 * l / lc);
   }
 }

 std::vector<Healpix_Map<double> > transformBspline_hp(Healpix_Map<double>& map, int m_nbScales) {
   int m_npix = map.Npix();
   std::vector<Healpix_Map<double> > band;
   band.push_back(map);
   for (int b=0; b<m_nbScales-1; b++) {
     Healpix_Map<double> map_out = smoothBspline_hp(map, b);
     band.push_back(map_out);
     for (int it = 0; it<m_npix; it++) {
       band[b][it] = band[b][it] - band[b+1][it];
       //if (isnan(band[b][it])) {
        // band[b][it]=band[b][it];
       //}
     }
   }
   return band;
 }

 Healpix_Map<double> smoothBspline_hp(Healpix_Map<double>& map, int scale) {
   int m_nside = map.Nside();
   int m_nlmax = 3*m_nside-1;
   Alm<xcomplex<double> > map_lm(m_nlmax, m_nlmax);
   map_lm.SetToZero();
   arr<double> weight;
   weight.alloc(2*m_nside);
   weight.fill(1.);

   map2alm_iter (map, map_lm, 3, weight);
   double lc = m_nlmax*pow(0.5, scale);
   //int nalm = int(m_nlmax*(m_nlmax+1)/2+(m_nlmax+1));

   double* Filter = new double[m_nlmax+1];
   std::fill_n(Filter, m_nlmax+1, 0);
   //double* Filter = new double[nalm];
   //std::fill_n(Filter, nalm, 0);
   getFilter(Filter, lc, m_nlmax);

   Alm<xcomplex<double> > smoothmap_lm(m_nlmax, m_nlmax);

   smoothmap_lm.SetToZero();
   for (int l=0; l<=m_nlmax; l++) {
     for (int m=0; m<=l; m++) {
       //smoothmap_lm(l,m) = Filter[getidx_lm(m_nlmax, l, m)] * map_lm(l,m);
       smoothmap_lm(l,m) = Filter[l] * map_lm(l,m);
     }
   }
   Healpix_Map<double> smooth_map;
   smooth_map.SetNside(m_nside, RING);
   smooth_map.fill(0.);
   alm2map(smoothmap_lm, smooth_map);
   delete[] Filter;
   weight.dealloc();
   return smooth_map;
 }

 Healpix_Map<double> reconsBspline_hp(std::vector<Healpix_Map<double> >& band) {
   int m_nside = band[0].Nside();
   //int m_nlmax = 3*m_nside-1;
   int m_npix = band[0].Npix();
   Healpix_Map<double> map_out;
   map_out.SetNside(m_nside, RING);
   map_out.fill(0.);
   for (size_t b = 0; b<band.size(); b++) {
     for (int it = 0; it<m_npix; it++) {
       map_out[it] = map_out[it] + band[b][it];
     }
   }
   return map_out;
 }

   } /* namespace Spherical */
  } /* namespace TwoDMass */
 } /* namespace WeakLensing */
} /* namespace Euclid */
