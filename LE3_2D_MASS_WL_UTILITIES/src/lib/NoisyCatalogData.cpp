/**
 * @file src/lib/NoisyCatalogData.cpp
 * @date 10/13/20
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

#include "LE3_2D_MASS_WL_UTILITIES/NoisyCatalogData.h"
#include <cmath>

static Elements::Logging logger = Elements::Logging::getLogger("RandomizeData");

namespace Euclid {
 namespace WeakLensing {
  namespace TwoDMass {

  NoisyCatalogData::NoisyCatalogData(){}

  std::vector<std::vector<double> > NoisyCatalogData::create_noisy_data(std::vector<std::vector<double> > inputData){
   std::vector<std::vector<double> > Output;
   std::vector <double> gamma1;
   std::vector <double> gamma2;
   std::vector <double> gamma1_n;
   std::vector <double> gamma2_n;
   unsigned int galCount = 0;
   gamma1 = inputData[3];
   gamma2 = inputData[4];
   galCount = inputData[0].size();
   //logger.info()<<"Total Galaxy count: "<<galCount;

   for (unsigned int i = 0; i< galCount; i++){
    float theta, radius, temp;
    temp = pow(gamma1[i], 2) + pow(gamma2[i], 2);
    radius = sqrt(temp);
    theta = M_PI * ((double) rand() / (RAND_MAX));
    gamma1_n.push_back(radius*cos(2*theta));
    gamma2_n.push_back(radius*sin(2*theta));
   }
   Output.push_back(inputData[0]);
   Output.push_back(inputData[1]);
   Output.push_back(inputData[2]);
   Output.push_back(gamma1_n);
   Output.push_back(gamma2_n);
   Output.push_back(inputData[5]);
   Output.push_back(inputData[6]);
   return Output;
  }

  } /* namespace TwoDMass */
 } /* namespace WeakLensing */
} /* namespace Euclid */
