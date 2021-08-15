/**
 * @file src/lib/MapMaker.cpp
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

#include "LE3_2D_MASS_WL_CARTESIAN/MapMaker.h"
#include "LE3_2D_MASS_WL_CARTESIAN/Projection.h"

#include <cmath>
using namespace Euclid::WeakLensing::TwoDMass;
using LE3_2D_MASS_WL_CARTESIAN::Projection;
static Elements::Logging logger = Elements::Logging::getLogger("MapMaker");
namespace LE3_2D_MASS_WL_CARTESIAN {

MapMaker::MapMaker(std::vector<std::vector<double> >& inData,
          LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam):
          inputData(inData), cartesianParam(cartesianParam)
{ }

ShearMap* MapMaker::getShearMap(LE3_2D_MASS_WL_CARTESIAN::CoordinateBound& CB){
 std::pair<long, double*> arrayPair =
           getMap(Euclid::WeakLensing::TwoDMass::mapType::shearMap, CB);
  // If there is no output return a nullptr
  if (arrayPair.second==nullptr) {
    return nullptr;
  }
  if (arrayPair.first==0) {
    delete [] arrayPair.second;
    arrayPair.second = nullptr;
    return nullptr;
  }
 int xbin = cartesianParam.getXaxis();
 int ybin = xbin;
  // Otherwise create the map
  ShearMap *myShearMap = new ShearMap(arrayPair.second, xbin, ybin, 3, CB, arrayPair.first);

  // Destroy the memory of the array
  delete [] arrayPair.second;
  arrayPair.second = nullptr;

  // And finally return the map
  return myShearMap;
}

ConvergenceMap* MapMaker::getConvMap(LE3_2D_MASS_WL_CARTESIAN::CoordinateBound& CB){
 std::pair<long, double*> arrayPair =
          getMap(Euclid::WeakLensing::TwoDMass::mapType::convMap, CB);
  // If there is no output return a nullptr
  if (arrayPair.second==nullptr) {
    return nullptr;
  }
  if (arrayPair.first==0) {
    delete [] arrayPair.second;
    arrayPair.second = nullptr;
    return nullptr;
  }
 int xbin = cartesianParam.getXaxis();
 int ybin = xbin;
  // Otherwise create the map
  ConvergenceMap *myConvMap = new ConvergenceMap(arrayPair.second, xbin, ybin, 3, CB, arrayPair.first);

  // Destroy the memory of the array
  delete [] arrayPair.second;
  arrayPair.second = nullptr;

  // And finally return the map
  return myConvMap;
}

std::pair<long, double*> MapMaker::getMap(const Euclid::WeakLensing::TwoDMass::mapType type,
         LE3_2D_MASS_WL_CARTESIAN::CoordinateBound& CB){
 //int xbin, ybin;
 int xbin =  cartesianParam.getXaxis();
 int ybin = cartesianParam.getYaxis();
 LE3_2D_MASS_WL_CARTESIAN::Projection Projection;
 double binXSize, binYSize;
 //double raMin, raMax, decMin, decMax;

 double raMin = CB.getRaMin();
 double raMax = CB.getRaMax();
 double decMin = CB.getDecMin();
 double decMax = CB.getDecMax();
 //logger.info()<<"ra min and max: "<<raMin<<" "<<raMax;
 //logger.info()<<"dec min and max: "<<decMin<<" "<<decMax;
 //logger.info()<<"z min and max: "<<CB.getZMin()<<" "<<CB.getZMax();
 //logger.info()<<"Param z min and max: "<<cartesianParam.getZMin()<<" "<<cartesianParam.getZMax();
 double ra0 = 0.5*(raMax + raMin);
 double dec0 = 0.5*(decMax + decMin);
 double raRange = (raMax - raMin);
 double decRange = (decMax - decMin);

 std::pair<double, double> xyMin = Projection.getGnomonicProjection(raMin, decMin, ra0, dec0);
 //xMin = xyMin.first; yMin = xyMin.second
 // In case a square map is not expected, perform some basic radec projection ranges selection
 if (cartesianParam.getSquareMap() == false) {
  std::pair<double, double> xyMax = Projection.getGnomonicProjection(raMax, decMax, ra0, dec0);
 //xMax = xyMax.first; yMin = xyMin.second
  //logger.info()<<"xMin:"<<xyMin.first<<std::endl;
  //logger.info()<<"xMax:"<<xyMax.first<<std::endl;
  // Define the bins sizes
  binXSize = (xyMax.first - xyMin.first)/xbin;
  binYSize = (xyMax.second - xyMin.second)/ybin;
 } else {
   std::pair<double, double> raDec1 =
          Projection.getInverseGnomonicProjection(-0.5*raRange*M_PI/180., -0.5*decRange*M_PI/180., ra0, dec0);
   std::pair<double, double> raDec2 =
          Projection.getInverseGnomonicProjection(0, -0.5*decRange*M_PI/180., ra0, dec0);
   std::pair<double, double> raDec3 =
          Projection.getInverseGnomonicProjection(0.5*raRange*M_PI/180., -0.5*decRange*M_PI/180., ra0, dec0);
   std::pair<double, double> raDec4 =
          Projection.getInverseGnomonicProjection(0.5*raRange*M_PI/180., 0.5*decRange*M_PI/180., ra0, dec0);
   std::pair<double, double> raDec5 =
          Projection.getInverseGnomonicProjection(0, 0.5*decRange*M_PI/180., ra0, dec0);
   std::pair<double, double> raDec6 =
          Projection.getInverseGnomonicProjection(-0.5*raRange*M_PI/180., 0.5*decRange*M_PI/180., ra0, dec0);

   // Get the min and max values of ra and dec according to the geometrical effects of projection
   raMin = raDec1.first < raDec6.first ? raDec1.first : raDec6.first;
   decMin = raDec1.second < raDec2.second ? raDec1.second : raDec2.second;
   raMax = raDec3.first > raDec4.first ? raDec3.first : raDec4.first;
   decMax = raDec4.second > raDec5.second ? raDec4.second : raDec5.second;

/*      logger.info()<<"ra dec 1: "<<raDec1.first<<" "<<raDec1.second;
      logger.info()<<"ra dec 2: "<<raDec2.first<<" "<<raDec2.second;
      logger.info()<<"ra dec 3: "<<raDec3.first<<" "<<raDec3.second;
      logger.info()<<"ra dec 4: "<<raDec4.first<<" "<<raDec4.second;
      logger.info()<<"ra dec 5: "<<raDec5.first<<" "<<raDec5.second;
      logger.info()<<"ra dec 6: "<<raDec6.first<<" "<<raDec6.second;
      logger.info()<<"ra min and max: "<<raMin<<" "<<raMax;
      logger.info()<<"dec min and max: "<<decMin<<" "<<decMax;
*/
      // Define the bins sizes
      binXSize = raRange*M_PI/180./xbin;
      binYSize = decRange*M_PI/180./ybin;
 //logger.info()<<"binXSize: "<<binXSize<<" 	binYSize: "<<binYSize;
 }
 //double *Array = initializeArray<double>(xbin, ybin, 2);
 //double *countArray = initializeArray<double>(xbin, ybin, 1);
 double *Array = new double[xbin*ybin*2];
 std::fill_n(Array, xbin*ybin*2, 0);
 double *countArray = new double[xbin*ybin*1];
 std::fill_n(countArray, xbin*ybin*1, 0);
 unsigned int galCount = 0, selGalCount = 0; //count for galaxies watched and selected

 for (unsigned int i=0; i<inputData[0].size(); i++){
  //std::vector<double> weight;
  double weight = 1.0; //if there is weight column in catalogue then use that instead of this
  if (inputData[6].empty()==false){
   weight = inputData[6][i];
  }
  if (inputData[1][i]>= decMin && inputData[1][i]<= decMax){
   if (inputData[0][i]>= raMin && inputData[0][i]<= raMax){
    if (inputData[5][i]>= CB.getZMin() && inputData[5][i]<= CB.getZMax()) {
     int tmpx, tmpy;
     if (cartesianParam.getSquareMap() == true) {
      // project the selected radec on gnomonic plan
      std::pair<double, double> tmpXY = Projection.getGnomonicProjection(inputData[0][i], inputData[1][i], ra0, dec0);
      //  calculate where it is on the binning
      tmpx = int(floor((tmpXY.first+0.5*raRange*M_PI/180.)/binXSize));
      tmpy = int(floor((tmpXY.second+0.5*decRange*M_PI/180.)/binYSize));
     } else {
      std::pair<double, double> tmpXY = Projection.getGnomonicProjection(inputData[0][i], inputData[1][i], ra0, dec0);
      tmpx = int(floor((tmpXY.first-xyMin.first)/binXSize));
      tmpy = int(floor((tmpXY.second-xyMin.second)/binYSize));
     }
    // logger.info()<< "temp: "<<tmpx;
     if (tmpx>=0 && tmpx<int(xbin) && tmpy>=0 && tmpy<int(ybin)){
      if (type == Euclid::WeakLensing::TwoDMass::mapType::shearMap){
       // Apply correction for the projection
       std::pair<double, double> xy1 = Projection.getGnomonicProjection(inputData[0][i], inputData[1][i], ra0, dec0);
       std::pair<double, double> xy2 =
                                Projection.getGnomonicProjection(inputData[0][i], (inputData[1][i]+0.01), ra0, dec0);
       double rotationAngle = -atan((xy2.first-xy1.first)/(xy2.second-xy1.second));
       double gamma1cor = -(inputData[3][i]*cos(2*rotationAngle)-inputData[4][i]*sin(2*rotationAngle));
       double gamma2cor = -(inputData[3][i]*sin(2*rotationAngle)+inputData[4][i]*cos(2*rotationAngle));
       Array[tmpy*xbin + tmpx] += gamma1cor*weight;
       Array[xbin*ybin + tmpy*xbin + tmpx] += gamma2cor*weight;
      }
      if (type == Euclid::WeakLensing::TwoDMass::mapType::convMap){
       Array[tmpy*xbin + tmpx] += inputData[2][i]*weight;
      }
      countArray[tmpy*xbin + tmpx] += weight;
      selGalCount+=weight;
   //  }
   // }
  // }
  } } } }
   galCount+=weight;
 }
 logger.info()<<"number of galaxies selected: "<<selGalCount;
 logger.info()<<"over the total number of galaxies: "<<galCount;
 // Normalize the values to have the mean shear in each bin
 for (int i = 0; i<xbin*ybin; i++) {
  if (countArray[i]>1) {
   Array[i] /= countArray[i];
   Array[i + xbin*ybin] /= countArray[i];
  //}
 } }
 if (Array == nullptr){
   logger.info()<< "INFO: data array is empty. Unable to create Shear map";
   //return false;
   exit (EXIT_FAILURE);
 }

 if (selGalCount == 0){
   logger.info()<< "INFO: no galaxies are selected";
   //return false;
   exit (EXIT_FAILURE);
 }
 //double *outArray = initializeArray<double>(xbin, ybin, 3);
 double *outArray = new double[xbin*ybin*3];
 std::fill_n(outArray, xbin*ybin*3, 0);
 //int sizeArray = xbin * ybin * 2;
 //int sizeCountArray = xbin * ybin *1;
 std::copy(Array, Array + (xbin * ybin * 2), outArray);
 std::copy(countArray, countArray + (xbin * ybin *1), outArray + (xbin * ybin * 2));
 // Delete dynamic arrays
 delete [] countArray;
 delete [] Array;
 return std::pair<long, double*> (selGalCount, outArray);

}

}  // namespace LE3_2D_MASS_WL_CARTESIAN
