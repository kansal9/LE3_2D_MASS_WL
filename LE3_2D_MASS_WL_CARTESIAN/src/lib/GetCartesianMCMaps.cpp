/**
 * @file src/lib/GetCartesianMCMaps.cpp
 * @date 10/13/20
 * @author vkansal
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

#include "LE3_2D_MASS_WL_CARTESIAN/GetCartesianMCMaps.h"
#include "LE3_2D_MASS_WL_CARTESIAN/MapMaker.h"
#include "LE3_2D_MASS_WL_CARTESIAN/MassMapping.h"
#include "LE3_2D_MASS_WL_UTILITIES/NoisyCatalogData.h"
#include "LE3_2D_MASS_WL_CARTESIAN/Matrix.h"
#include "ElementsKernel/Temporary.h"
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

using namespace Euclid::WeakLensing::TwoDMass;
 // handle on created path names
 boost::filesystem::path temp_path;
static Elements::Logging logger = Elements::Logging::getLogger("CartesianMCMaps");

namespace LE3_2D_MASS_WL_CARTESIAN {

 GetCartesianMCMaps::GetCartesianMCMaps (std::vector<std::vector<double> >& inputData,
      LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam, LE3_2D_MASS_WL_CARTESIAN::CoordinateBound& CB):
             m_inData(inputData), m_cartesianParam(cartesianParam), m_CB(CB)
{ }

 LE3_2D_MASS_WL_CARTESIAN::ShearMap* GetCartesianMCMaps::getDeNoisedShearMap() {
  LE3_2D_MASS_WL_CARTESIAN::ShearMap *DeNoisedShearMap;

  Elements::TempDir one;
  boost::filesystem::path temp_path = one.path();
 // Save this map back in a new file$
  boost::filesystem::path shearMap = temp_path / "ShearMapFits.fits";

 // Object to perform extraction
  Euclid::WeakLensing::TwoDMass::CartesianKS::CartesianAlgoKS CartesainAlgo(m_cartesianParam);
  CartesainAlgo.extractShearMap(shearMap.native(), m_inData, m_CB);

  logger.info() << "size cat: " << m_inData[1].size();

  boost::filesystem::path convMap = temp_path / "convMapFits.fits";
//  performKSMassMapping(shearMap_MC.native(), convMapFits_MC.native());
  LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap *convMapFits;
  MassMapping mass(m_cartesianParam);

  convMapFits = mass.getSheartoConv(shearMap.native());

  if (fabs(m_cartesianParam.getSigmaGauss())>0.001){//Apply gaussian filter on the map
   convMapFits->applyGaussianFilter(m_cartesianParam.getSigmaGauss()); //this will update contents by applying filter
  }

  // Writing Convergence Map
  if (convMapFits != nullptr) {
   convMapFits->writeMap(convMap.native(), m_cartesianParam);
  }

  // Perform inverse mass mapping
  MassMapping inversion( m_cartesianParam);
  DeNoisedShearMap = inversion.getConvtoShear(convMap.native());

  DeNoisedShearMap->computeReducedShear(*convMapFits);  

  delete convMapFits;
  convMapFits = nullptr;

  return DeNoisedShearMap;
 }

 std::vector<LE3_2D_MASS_WL_CARTESIAN::ShearMap*> GetCartesianMCMaps::getNoisedShearMaps() {
  int Niter = m_cartesianParam.getNSamples();
  if (Niter <= 0) {
   Niter = 1;
  }
  std::vector<LE3_2D_MASS_WL_CARTESIAN::ShearMap*> ShearMapList;

  for (int i = 0; i !=Niter; ++i) {
    std::vector<std::vector<double> > RanData;
    NoisyCatalogData randomise;
    RanData = randomise.create_noisy_data(m_inData);
    MapMaker mapMaker(RanData, m_cartesianParam);
    LE3_2D_MASS_WL_CARTESIAN::ShearMap *NoisedShearMap;
    NoisedShearMap = mapMaker.getShearMap(m_CB);
    ShearMapList.push_back(NoisedShearMap);
   // const std::string map = "noisedShearMap_" + std::to_string(i) + ".fits";
   // NoisedShearMap->writeMap(map, m_cartesianParam);
  }
   return ShearMapList;
 }

  bool GetCartesianMCMaps::performAddition(LE3_2D_MASS_WL_CARTESIAN::ShearMap& DenoisedShearMap,
                        LE3_2D_MASS_WL_CARTESIAN::ShearMap& NoisedShearMap, const std::string& ShearMap ) {

   LE3_2D_MASS_WL_CARTESIAN::Matrix ShearE(m_cartesianParam.getXaxis(), m_cartesianParam.getYaxis());
   LE3_2D_MASS_WL_CARTESIAN::Matrix ShearB(m_cartesianParam.getXaxis(), m_cartesianParam.getYaxis());

    for (int j=0; j<m_cartesianParam.getYaxis(); j++) {
     for (int i=0; i<m_cartesianParam.getXaxis(); i++) {
       double DenoisedE = DenoisedShearMap.getBinValue(i, j, 0);
       double DenoisedB = DenoisedShearMap.getBinValue(i, j, 1);
       double NoisedE = NoisedShearMap.getBinValue(i, j, 0);
       double NoisedB = NoisedShearMap.getBinValue(i, j, 1);
       double r1 = (1. + ((DenoisedE * NoisedE) - ((-1. * DenoisedB) * NoisedB)));
       double i1 = (DenoisedE * NoisedB) + ((-1. * DenoisedB) * NoisedE);
       double r2 = NoisedE + DenoisedE;
       double i2 = NoisedB + DenoisedB;
       double valE = ((r2*r1) + (i2*i1))/(pow(r1, 2) + pow(i1, 2));
       double valB = ((i2*r1) - (i1*r2))/(pow(r1, 2) + pow(i1, 2));
       ShearE.setValue(i, j, valE);
       ShearB.setValue(i, j, valB);
     }
    }

   double *mapArray = new double[m_cartesianParam.getXaxis()*m_cartesianParam.getYaxis()*2];
   // Fill the map array
   for (int i=0; i<m_cartesianParam.getXaxis(); i++){
    for (int j=0; j<m_cartesianParam.getYaxis(); j++) {
     mapArray[j*m_cartesianParam.getXaxis() +i] = ShearE.getValue(i, j);
     if (m_cartesianParam.getForceBMode() == true){
     mapArray[m_cartesianParam.getXaxis()*m_cartesianParam.getYaxis() + j*m_cartesianParam.getXaxis() +i] = 0;
     } else {
     mapArray[m_cartesianParam.getXaxis()*m_cartesianParam.getYaxis() + j*m_cartesianParam.getXaxis() +i]
                                                            = ShearB.getValue(i, j);
     }
    }
   }
   LE3_2D_MASS_WL_CARTESIAN::ShearMap shearMap(mapArray, m_cartesianParam.getXaxis(), m_cartesianParam.getYaxis(), 2,
                                               NoisedShearMap.getNumberOfGalaxies());
   shearMap.writeMap(ShearMap, m_cartesianParam);
   delete [] mapArray;
   mapArray = nullptr;
   return true;
  }

}  // namespace LE3_2D_MASS_WL_CARTESIAN
