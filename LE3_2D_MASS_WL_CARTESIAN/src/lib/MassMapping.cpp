/**
 * @file src/lib/MassMapping.cpp
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

#include "LE3_2D_MASS_WL_CARTESIAN/MassMapping.h"
#include "LE3_2D_MASS_WL_CARTESIAN/InpaintingAlgo.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <pthread.h>
#include <typeinfo>
#include <cstring>
#include <climits>
#include <cmath>
#include <vector>
#include <list>
#include <map>
#include <fftw3.h>
#include "boost/multi_array.hpp"

using namespace Euclid::WeakLensing::TwoDMass;
using LE3_2D_MASS_WL_CARTESIAN::InpaintingAlgo;
static Elements::Logging logger = Elements::Logging::getLogger("MassMapping");
namespace LE3_2D_MASS_WL_CARTESIAN {

MassMapping::MassMapping(LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam): m_convMap(nullptr),
                         m_shearMap(nullptr), m_cartesianParam(cartesianParam) {
}

ConvergenceMap* MassMapping::getSheartoConv(const std::string& map){
  m_shearMap = new ShearMap(map);
  /*if (m_shearMap==nullptr) {
   return false;
  }*/
  if (fabs(m_cartesianParam.getSigmaGauss())>0.001){//Apply gaussian filter on the map
   m_shearMap->applyGaussianFilter(m_cartesianParam.getSigmaGauss()); //this will update contents by applying filter
  }
  /*if (m_cartesianParam.get_removeOffset()){
   m_shearMap->removeOffset(m_shearMap->getMeanValues());
   std::vector<double> tmpMean = m_shearMap->getMeanValues();
   logger.info()<<"mean value of G1: "<<tmpMean[0];
   logger.info()<<"mean value of G2: "<<tmpMean[1];
  }*/

  if (m_cartesianParam.get_addBorders()){
   m_shearMap->add_borders();
  }
  //std::pair<unsigned int, double*> outMap = Map->massMap();
  /* logger.info()<<"Nb Gal: "<< m_shearMap->getNumberOfGalaxies();
   logger.info()<<"RaMax: "<< m_shearMap->getCoordinateBound().getRaMax();
   logger.info()<<"X dim: "<< m_shearMap->getXdim();
   logger.info()<<"Y dim: "<< m_shearMap->getYdim();
   logger.info()<<"Z dim: "<< m_shearMap->getZdim();*/
   LE3_2D_MASS_WL_CARTESIAN::CoordinateBound mCB = m_shearMap->getCoordinateBound();
   ConvergenceMap convergence = m_shearMap->getConvMap();
   double *MapArr = new double[m_shearMap->getXdim()*m_shearMap->getYdim()*m_shearMap->getZdim()];
   convergence.getArray(MapArr);
   //m_convMap = new ConvergenceMap(m_shearMap->getConvMap());
   m_convMap = new ConvergenceMap(MapArr, m_shearMap->getXdim(), m_shearMap->getYdim(), m_shearMap->getZdim(),
                                  mCB, m_shearMap->getNumberOfGalaxies());
   delete [] MapArr;
   MapArr = nullptr;
   delete m_shearMap;
   m_shearMap = nullptr;
   return m_convMap;
}

ShearMap* MassMapping::getConvtoShear(const std::string& map){
   m_convMap = new ConvergenceMap(map);
   if (fabs(m_cartesianParam.getSigmaGauss())>0.001){//Apply gaussian filter on the map
    m_convMap->applyGaussianFilter(m_cartesianParam.getSigmaGauss()); //this will update contents by applying filter
   }
   /*if (m_cartesianParam.get_removeOffset() == 1){
    m_convMap->removeOffset(m_convMap->getMeanValues());
    std::vector<double> tmpMean = m_convMap->getMeanValues();
    logger.info()<<"mean value of G1: "<<tmpMean[0];
    logger.info()<<"mean value of G2: "<<tmpMean[1];
   }*/
   //m_shearMap = new ShearMap(m_convMap->getShearMap());
   LE3_2D_MASS_WL_CARTESIAN::CoordinateBound mCB = m_convMap->getCoordinateBound();
   ShearMap shear = m_convMap->getShearMap();
   double *MapArr = new double[m_convMap->getXdim()*m_convMap->getYdim()*m_convMap->getZdim()];
   shear.getArray(MapArr);
   m_shearMap = new ShearMap(MapArr, m_convMap->getXdim(), m_convMap->getYdim(), m_convMap->getZdim(),
                                  mCB, m_convMap->getNumberOfGalaxies());
   delete [] MapArr;
   MapArr = nullptr;
   delete m_convMap;
   m_convMap = nullptr;
   return m_shearMap;
}
}  // namespace LE3_2D_MASS_WL_CARTESIAN
