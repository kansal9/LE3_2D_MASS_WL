/**
 * @file src/lib/MatrixProcess.cpp
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

#include "LE3_2D_MASS_WL_CARTESIAN/MatrixProcess.h"
#include "fftw3.h"
#include "math.h"
#include <iostream>
static Elements::Logging logger = Elements::Logging::getLogger("MatrixProcess");
namespace LE3_2D_MASS_WL_CARTESIAN {

MatrixProcess::MatrixProcess(unsigned int sizeXaxis, unsigned int sizeYaxis):
                            m_sizeXaxis(sizeXaxis), m_sizeYaxis(sizeYaxis) { }

Matrix MatrixProcess::performDCT(Matrix input) {
 // Create an output image
 Matrix DCToutput(m_sizeXaxis, m_sizeYaxis);
 // Create the plan on which to perform the transform
 fftw_plan DCTplan;

 #pragma omp critical
  {
    DCTplan = fftw_plan_r2r_2d(m_sizeXaxis, m_sizeYaxis, input.getArray(), DCToutput.getArray(),
                               FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
  }

 // Perform the transformation
 fftw_execute(DCTplan);

 // Free memory
 fftw_destroy_plan(DCTplan);

 // Rescale the output
 double dctFactor = 2*sqrt(m_sizeXaxis*m_sizeYaxis);
 return DCToutput.multiply(1./dctFactor);
}

Matrix MatrixProcess::performIDCT(Matrix input) {
 // Create an output image
 Matrix output(m_sizeXaxis, m_sizeYaxis);
 // Create the plan on which to perform the inverse transform
 fftw_plan IDCTplan;

 #pragma omp critical
 {
 IDCTplan = fftw_plan_r2r_2d(m_sizeXaxis, m_sizeYaxis, input.getArray(), output.getArray(),
                                FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
 }
 // Perform the transformation
 fftw_execute(IDCTplan);
 // Free memory
 fftw_destroy_plan(IDCTplan);
 // Rescale the output
 double dctFactor = 2*sqrt(m_sizeXaxis*m_sizeYaxis);
 return output.multiply(1./dctFactor);
}

Matrix MatrixProcess::performDCT(Matrix input, unsigned int blockSizeX, unsigned int blockSizeY, bool forward) {
 // Create an output image
 Matrix output(m_sizeXaxis, m_sizeYaxis);
 MatrixProcess myNewMatrix(blockSizeX, blockSizeY);
 
 for (unsigned iblock=0; iblock<m_sizeXaxis/blockSizeX; iblock++) {
  for (unsigned jblock=0; jblock<m_sizeYaxis/blockSizeY; jblock++) {
   // Fill the block map
   Matrix blockInput(blockSizeX, blockSizeY);
   for (size_t i=0; i<blockSizeX; i++) {
    for (size_t j=0; j<blockSizeY; j++) {
     blockInput.setValue(i, j, input.getValue(iblock*blockSizeX+i, jblock*blockSizeY+j));
    }
   }
   // Perform DCT on that block map
   Matrix blockOutput = forward ? myNewMatrix.performDCT(blockInput) : myNewMatrix.performIDCT(blockInput);
   // Fill the global output with that block map
   for (size_t i=0; i<blockSizeX; i++) {
    for (size_t j=0; j<blockSizeY; j++) {
     output.setValue(iblock*blockSizeX+i, jblock*blockSizeY+j, blockOutput.getValue(i, j));
    }
   }
  }
 }
 // return the global output
 return output;
}

std::vector<Matrix> MatrixProcess::transformBspline(Matrix input, unsigned int nbScales) {
 // Create a vector of images
 std::vector<Matrix> band;
 // Add the first image to this vector
 band.push_back(input);
 // Loop over the number of scales
 for (size_t step=0; step<nbScales-1; step++) {
  // Apply the b spline transfo with algorithm a trous
  Matrix imageOut = smoothBspline(band[step], step);
  // Add the output image to the vector
  band.push_back(imageOut);
  band[step] = band[step].substract(band[step+1]);
 }
 return band;
}

Matrix MatrixProcess::smoothBspline(Matrix input, unsigned int stepTrou) {
 // Define some values for the transform
 float h0 = 3./8.;
 float h1 = 1./4.;
 float h2 = 1./16.;
 int step = int(pow(2., double(stepTrou))+0.5);
 // Create an intermediate image
 Matrix tmpImage(m_sizeXaxis, m_sizeYaxis);

 for (size_t i=0; i<m_sizeXaxis; i++) {
  for (size_t j=0; j<m_sizeYaxis; j++) {
   tmpImage.setValue(i, j, h0*input.getValue(i, j)
                             +h1*(input.getValue(i, int(j-step))+input.getValue(i, j+step))
                             +h2*(input.getValue(i, int(j-2*step))+input.getValue(i, j+2*step)));
  }
 }

 // Create the output image
 Matrix imageOut(m_sizeXaxis, m_sizeYaxis);
 for (size_t i=0; i<m_sizeXaxis; i++) {
  for (size_t j=0; j<m_sizeYaxis; j++) {
   imageOut.setValue(i, j, h0*tmpImage.getValue(i, j)
                             +h1*(tmpImage.getValue(int(i-step), j)+tmpImage.getValue(i+step, j))
                             +h2*(tmpImage.getValue(int(i-2*step), j)+tmpImage.getValue(i+2*step, j)));
  }
 }
 return imageOut;
}

Matrix MatrixProcess::reconsBspline(std::vector<Matrix> band) {
 // Create the output image
  Matrix imageOut(band[0].getXdim(), band[0].getYdim());

 // Just add the images from all the scales
 for (size_t i=0; i<band.size(); i++) {
  imageOut = imageOut.add(band[i]);
 }
 return imageOut;
}
}  // namespace LE3_2D_MASS_WL_CARTESIAN
