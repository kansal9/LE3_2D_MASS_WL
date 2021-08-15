/**
 * @file src/lib/ConvergenceMap.cpp
 * @date 05/19/20
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

#include "LE3_2D_MASS_WL_CARTESIAN/ConvergenceMap.h"

using namespace Euclid::WeakLensing::TwoDMass;
static Elements::Logging logger = Elements::Logging::getLogger("ConvergenceMap");
namespace LE3_2D_MASS_WL_CARTESIAN {
 ConvergenceMap::ConvergenceMap(double* array, LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam,
                                int nbGalaxies): GetMap(array, cartesianParam, nbGalaxies)
{
}

 ConvergenceMap::ConvergenceMap(double* array, int sizeXaxis, int sizeYaxis,  int sizeZaxis,
                                int nbGalaxies): GetMap(array, sizeXaxis, sizeYaxis, sizeZaxis, nbGalaxies)
{
}
 ConvergenceMap::ConvergenceMap(double* array, int sizeXaxis, int sizeYaxis,  int sizeZaxis,
             LE3_2D_MASS_WL_CARTESIAN::CoordinateBound& coordBound, int nbGalaxies):
             GetMap(array, sizeXaxis, sizeYaxis, sizeZaxis, coordBound, nbGalaxies)
{
}
 ConvergenceMap::ConvergenceMap(const std::string& filename): GetMap(filename){}

 ConvergenceMap::ConvergenceMap(GetMap const& copyMap):GetMap(copyMap) {}


ShearMap ConvergenceMap::getShearMap(){
  double fftFactor = 1.0/sizeXaxis/sizeYaxis;

  // Create the complex maps
  fftw_complex* kappa_complex = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*sizeXaxis*sizeYaxis);
  fftw_complex* fft_kappa_complex = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*sizeXaxis*sizeYaxis);
  fftw_complex* Psi_complex = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *sizeXaxis*sizeYaxis);
  fftw_complex* gamma_complex = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *sizeXaxis*sizeYaxis);
  fftw_complex* fft_gamma_complex = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *sizeXaxis*sizeYaxis);

  // Create the plans for transformation
  fftw_plan plan_k_forward;
  fftw_plan plan_g_backward;

  #pragma omp critical
  {
    plan_k_forward = fftw_plan_dft_2d(sizeXaxis, sizeYaxis, kappa_complex, fft_kappa_complex,
                                      FFTW_FORWARD,  FFTW_MEASURE);
    plan_g_backward = fftw_plan_dft_2d(sizeXaxis, sizeYaxis, fft_gamma_complex, gamma_complex,
                                       FFTW_BACKWARD,  FFTW_MEASURE);
  }

  // Fill the complex kappa map with convergence map values
  for ( int i=0; i<sizeXaxis; i++) {
    for (int j=0; j<sizeYaxis; j++) {
      kappa_complex[j*sizeXaxis +i][0] = (*m_mapValues)[i][j][0];
      kappa_complex[j*sizeXaxis +i][1] = (*m_mapValues)[i][j][1];
    }
  }

  // Perform the fourier transform of the complex convergence map
  fftw_execute(plan_k_forward);

  // Create the P factor
  for ( int i=0; i<sizeXaxis; i++) {
    for ( int j=0; j<sizeYaxis; j++) {
//        if (i+j == 0) continue;
      int ll1 = (double(i) <= double(sizeXaxis)/2. ? i : i - sizeXaxis);
      int ll2 = (double(j) <= double(sizeYaxis)/2. ? j : j - sizeYaxis);
      Psi_complex[j*sizeXaxis +i][0] = double(ll1*ll1-ll2*ll2)/double(ll1*ll1+ll2*ll2);
      Psi_complex[j*sizeXaxis +i][1] = double(2.*(ll1*ll2))/double(ll1*ll1+ll2*ll2);
   }
  }
  Psi_complex[0][0] = 0.0;
  Psi_complex[0][1] = 0.0;

  // Create the complex shear map in Fourier space
  for ( int i=0; i<sizeXaxis; i++) {
    for ( int j=0; j<sizeYaxis; j++) {
      fft_gamma_complex[j*sizeXaxis+i][0] = Psi_complex[j*sizeXaxis+i][0]*fft_kappa_complex[j*sizeXaxis+i][0]
                                             -Psi_complex[j*sizeXaxis+i][1]*fft_kappa_complex[j*sizeXaxis+i][1];
      fft_gamma_complex[j*sizeXaxis+i][1] = Psi_complex[j*sizeXaxis+i][0]*fft_kappa_complex[j*sizeXaxis+i][1]
                                             +Psi_complex[j*sizeXaxis+i][1]*fft_kappa_complex[j*sizeXaxis+i][0];
    }
  }

  // Perform inverse Fourier transform to get the shear map
  fftw_execute(plan_g_backward);

  // Fill the shear map
  double *gammaArray = new double[sizeXaxis*sizeYaxis*sizeZaxis];
  getArray(gammaArray);
  for ( int i=0; i<sizeXaxis; i++) {
    for ( int j=0; j<sizeYaxis; j++) {
      gammaArray[j*sizeXaxis +i] = gamma_complex[j*sizeXaxis +i][0]*fftFactor;
      gammaArray[sizeXaxis*sizeYaxis + j*sizeXaxis +i] = gamma_complex[j*sizeXaxis +i][1]*fftFactor;
   }
  }

  ShearMap gammaMap(gammaArray, sizeXaxis, sizeYaxis, sizeZaxis, nGalaxies);

  // free memory
  delete [] gammaArray;
  gammaArray = nullptr;

  fftw_destroy_plan(plan_k_forward);
  fftw_destroy_plan(plan_g_backward);

  fftw_free(gamma_complex);
  fftw_free(Psi_complex);
  fftw_free(kappa_complex);
  fftw_free(fft_gamma_complex);
  fftw_free(fft_kappa_complex);
  fftw_cleanup();
  return gammaMap;
}

}  // namespace LE3_2D_MASS_WL_CARTESIAN
