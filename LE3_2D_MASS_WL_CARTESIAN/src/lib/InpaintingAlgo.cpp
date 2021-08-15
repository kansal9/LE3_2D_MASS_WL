/**
 * @file src/lib/InpaintingAlgo.cpp
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

#include "LE3_2D_MASS_WL_CARTESIAN/InpaintingAlgo.h"
#include "LE3_2D_MASS_WL_CARTESIAN/MassMapping.h"

#include <map>
#include <fftw3.h>
#include <cmath>

using namespace Euclid::WeakLensing::TwoDMass;
static Elements::Logging logger = Elements::Logging::getLogger("Inpainting");
namespace LE3_2D_MASS_WL_CARTESIAN {

InpaintingAlgo::~InpaintingAlgo(){
 delete m_maskValues;
 m_maskValues = nullptr;
}

InpaintingAlgo::InpaintingAlgo(ShearMap &shearMap, ConvergenceMap &convMap,
           LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam):m_cartesianParam(cartesianParam),
           shearMap(shearMap), convMap(convMap), m_MP(shearMap.getXdim(), shearMap.getYdim()) {
 m_minThreshold = 0.0;
 m_maxThreshold = 0.5;
 typedef boost::multi_array<double, 3>::index index;
 // Retrieve the size of the maps
 Xaxis = shearMap.getXdim();
 Yaxis = shearMap.getYdim();
 Zaxis = shearMap.getZdim();
 nbScales = m_cartesianParam.getnbScales();
 if (nbScales==0) {
  nbScales = int(log(Xaxis)/log(2.))-2.;
 }
 logger.info()<<"axis dim: "<<Xaxis<<" "<<Yaxis<<" "<<Zaxis;
 logger.info()<<"number of scales: "<< nbScales;
 // Declare the mask
 m_maskValues = new boost::multi_array<int, 3>(boost::extents[Xaxis][Yaxis][Zaxis]);

 unsigned int count(0);

 // Initialize the mask
 for (index j = 0; j != Yaxis; ++j) {
  for (index i = 0; i != Xaxis; ++i) {
   if (fabs(shearMap.getBinValue(i, j, 0))<0.0000000001 && fabs(shearMap.getBinValue(i, j, 1))<0.0000000001)
   {
     (*m_maskValues)[i][j][0] = 0;
     (*m_maskValues)[i][j][1] = 0;
     count++;
   } else {
     (*m_maskValues)[i][j][0] = 1;
     (*m_maskValues)[i][j][1] = 1;
   }
  }
 }
 logger.info()<<"number of zeros: "<<count;
 logger.info()<<"over the number of pixels: "<<Xaxis*Yaxis*Zaxis;
}

InpaintingAlgo::InpaintingAlgo(const InpaintingAlgo& copy):
                  shearMap(copy.shearMap), convMap(copy.convMap),
                  m_MP(copy.m_MP), Xaxis(copy.Xaxis), Yaxis(copy.Yaxis), Zaxis(copy.Zaxis) {
 m_minThreshold = copy.m_minThreshold;
 m_maxThreshold = copy.m_maxThreshold;
 typedef boost::multi_array<double, 3>::index index;
 // Declare the mask
 m_maskValues = new boost::multi_array<int, 3>(boost::extents[Xaxis][Yaxis][Zaxis]);
 // Initialize the mask
 for (index j = 0; j != Yaxis; ++j) {
  for (index i = 0; i != Xaxis; ++i) {
   (*m_maskValues)[i][j][0] = (*(copy.m_maskValues))[i][j][0];
   (*m_maskValues)[i][j][1] = (*(copy.m_maskValues))[i][j][1];
  }
 }
}

ConvergenceMap* InpaintingAlgo::performInPaintingAlgo() {
 unsigned int nbIter = m_cartesianParam.getNInpaint();
 bool sigmaBounds = m_cartesianParam.getEqualVarPerScale();
 bool bModeZeros = m_cartesianParam.getForceBMode();

 logger.info()<<"number of inpainting iteration: "<< nbIter;
 logger.info()<<"status of Equal variance per scale is: "<< sigmaBounds;
 logger.info()<<"status of ForceBMode: "<< bModeZeros;
 // Create a copy of the conv map
 ConvergenceMap *kappaMapIter = new ConvergenceMap(convMap);

 double maxThreshold(m_maxThreshold);
//  double maxThreshold(-10.);
 double minThreshold(m_minThreshold);
//  double minThreshold(0.);

// logger.info()<<"Minimum threshold: "<<minThreshold;
// logger.info()<<"Maximum threshold: "<<maxThreshold;

 for (size_t iter = 0; iter<nbIter; iter++){
  logger.info()<<"iteration "<<iter<<" beginning";
  // Allocate memory for kappa and DCT kappa
  Matrix kappaE (Xaxis, Yaxis);
  Matrix kappaB (Xaxis, Yaxis);

 logger.info()<<"Initializing kappa map";
 for (size_t j=0; j<Yaxis; j++){
  for (size_t i=0; i<Xaxis; i++){
   kappaE.setValue(i, j, kappaMapIter->getBinValue(i, j, 0));
   kappaB.setValue(i, j, kappaMapIter->getBinValue(i, j, 1));
  }
 }
 // Perform the DCT
 Matrix DCTkappaE = m_MP.performDCT(kappaE);
 Matrix DCTkappaB = m_MP.performDCT(kappaB);

 // Update the threshold value with the max value at first iteration
 if (iter==0 && maxThreshold<=0.){
   maxThreshold = DCTkappaE.getMax();
 }
 double lambda = minThreshold + (maxThreshold-minThreshold)*(erfc(2.8*iter/nbIter));
 //if (iter==nbIter-1) {
 if (lambda < minThreshold || iter==nbIter-1) {
  lambda = minThreshold;
 }
 logger.info()<<"threshold: "<<lambda;

 // Keep in memory the 0 values of kappa
 double DCTkappaE0 = DCTkappaE.getValue(0, 0);
 double DCTkappaB0 = DCTkappaB.getValue(0, 0);

 // Cut all values below the threshold value
 DCTkappaE.applyThreshold(lambda);
 DCTkappaB.applyThreshold(lambda);

 // Restore 0 values
 DCTkappaE.setValue(0, 0, DCTkappaE0);
 DCTkappaB.setValue(0, 0, DCTkappaB0);

 // Perform the IDCT
 kappaE = m_MP.performIDCT(DCTkappaE);
 kappaB = m_MP.performIDCT(DCTkappaB);

 // Apply sigma boundaries
 if (sigmaBounds){
 kappaE = applyBoundariesOnWavelets(kappaE);
 }

 // Perform the inversion and apply the mask to get the final convergence map
 delete kappaMapIter;
 kappaMapIter = new ConvergenceMap(performInversionMask(kappaE, kappaB, bModeZeros));
  logger.info()<<"end of iteration "<<iter;
 }
 return kappaMapIter;
}

 ConvergenceMap* InpaintingAlgo::performInPaintingAlgo(unsigned int blockSizeX, unsigned int blockSizeY) {
 size_t nbIter = m_cartesianParam.getNInpaint();
 bool sigmaBounds = m_cartesianParam.getEqualVarPerScale();
 bool bModeZeros = m_cartesianParam.getForceBMode();
 // Create a copy of the conv map
 ConvergenceMap *kappaMapIter = new ConvergenceMap(convMap);
 float threshold1(0);
 float threshold2(0);
 for (size_t iter = 0; iter<nbIter; iter++) {
  logger.info()<<"iteration "<<iter<<" beginning";
  // Allocate memory for kappa and DCT kappa
  Matrix kappa_E(Xaxis, Yaxis);
  Matrix kappa_B(Xaxis, Yaxis);

 logger.info()<<"Initializing kappa map";
 for (size_t j=0; j<Yaxis; j++) {
  for (size_t i=0; i<Xaxis; i++){
   kappa_E.setValue(i, j, kappaMapIter->getBinValue(i, j, 0));
   kappa_B.setValue(i, j, kappaMapIter->getBinValue(i, j, 1));
  }
 }

 // Perform the DCT
 bool forward = true;
 Matrix DCTkappaE = m_MP.performDCT(kappa_E, blockSizeX, blockSizeY, forward);
 Matrix DCTkappaB = m_MP.performDCT(kappa_B, blockSizeX, blockSizeY, forward);

 // Update the threshold value with the max value at first iteration
 if (iter==0){
  threshold1 = DCTkappaE.getMax();
 }
 float lambda = threshold2 + (threshold1-threshold2)*(erfc(2.8*iter/nbIter));
 logger.info()<<"threshold: "<<lambda;
 // Cut all values below the threshold value excluding 0 values of each block
 for (size_t i=0; i<Xaxis; i++) {
  for (size_t j=0; j<Yaxis; j++){
   bool isNotFirstBlockValue = i%blockSizeX!=0 || j%blockSizeY!=0;
   if (isNotFirstBlockValue && fabs(DCTkappaE.getValue(i, j)) < lambda){
    DCTkappaE.setValue(i, j, 0);
   }
   if (isNotFirstBlockValue && fabs(DCTkappaB.getValue(i, j)) < lambda){
    DCTkappaB.setValue(i, j, 0);
   }
  }
 }
 // Perform the IDCT
 forward = false;
 kappa_E = m_MP.performDCT(DCTkappaE, blockSizeX, blockSizeY, forward);
 kappa_B = m_MP.performDCT(DCTkappaB, blockSizeX, blockSizeY, forward);

 // Apply sigma boundaries
 if (sigmaBounds){
  //    applyBoundaries(kappaE);
  kappa_E = applyBoundariesOnWavelets(kappa_E);
 }
 // Perform the inversion and apply the mask to get the final convergence map
 delete kappaMapIter;
 kappaMapIter = new ConvergenceMap(performInversionMask(kappa_E, kappa_B, bModeZeros));
 logger.info()<<"end of iteration "<<iter;
 }
 return kappaMapIter;
}

Matrix InpaintingAlgo::applyBoundariesOnWavelets(Matrix input) {

 std::vector<Matrix> myBand = m_MP.transformBspline(input, nbScales);
 for (int kScale=0; kScale<nbScales-1; kScale++) {
   double maskSigma=0.0, imSigma=0.0, maskCount=0.0, imCount=0.0;
   getMaskImageSigma(myBand[kScale], maskSigma, imSigma, maskCount, imCount);

   //logger.info()<<"maskSigma "<<maskSigma;
   //logger.info()<<"imSigma "<<imSigma;
   multiplyWaveletCoeff (myBand[kScale], maskSigma, imSigma, maskCount, imCount);
 }
 return m_MP.reconsBspline(myBand);
}

void InpaintingAlgo::getMaskImageSigma (LE3_2D_MASS_WL_CARTESIAN::Matrix& Image, double& maskSigma, double& imSigma,
                                        double& maskCount, double& imCount){

  double maskMean = 0.;
  double maskSquareMean = 0.;
  double imMean = 0.;
  double imSquareMean = 0.;
  for (size_t i=0; i<Xaxis; i++){
   for (size_t j=0; j<Yaxis; j++){
    double tmp = Image.getValue(i, j);
    if ((*m_maskValues)[i][j][0]==0){
     maskMean += tmp;
     maskSquareMean += tmp*tmp;
     maskCount++;
    } else {
     imMean += tmp;
     imSquareMean += tmp*tmp;
     imCount++;
    }
   }
  }
 maskSigma = sqrt((maskSquareMean/maskCount) - ((maskMean/maskCount)*(maskMean/maskCount)));
 imSigma = sqrt((imSquareMean/imCount) - ((imMean/imCount)*(imMean/imCount)));
}

void InpaintingAlgo::multiplyWaveletCoeff (LE3_2D_MASS_WL_CARTESIAN::Matrix& Image, double& maskSigma,
         double& imSigma, double& maskCount, double& imCount) {

   //if (maskSigma > imSigma) {
   if (maskSigma > imSigma*(1 + sqrt (sqrt (2. / (maskCount+1))))) {
    for (size_t i=0; i<Xaxis; i++){
     for (size_t j=0; j<Yaxis; j++){
      if ((*m_maskValues)[i][j][0]==0){
       if ((maskCount > 9) && (imCount > 9) && (maskSigma > 0)){
        Image.setValue(i, j, Image.getValue(i, j)*imSigma/maskSigma);
       }
      }
     }
    }
   }
}

ConvergenceMap InpaintingAlgo::performInversionMask(Matrix kappaE, Matrix kappaB, bool bModeZeros){
 // Allocate memory for the kappa array
 double *kappaArray = new double[Xaxis*Yaxis*Zaxis];
 // Fill the convergence map array
 for (size_t i=0; i<Xaxis; i++){
  for (size_t j=0; j<Yaxis; j++) {
   kappaArray[j*Xaxis +i] = kappaE.getValue(i, j);
   if ((*m_maskValues)[i][j][0]==0 && (bModeZeros == true)){
    kappaArray[Xaxis*Yaxis + j*Xaxis +i] = 0;
   } else {
    kappaArray[Xaxis*Yaxis + j*Xaxis +i] = kappaB.getValue(i, j);
   }
  }
 }

  ConvergenceMap kappaMap(kappaArray, Xaxis, Yaxis, Zaxis);

  // Free the memory
  delete [] kappaArray;
  kappaArray = nullptr;

  // Get the shear map from this convergence map
  ShearMap gammaMap = kappaMap.getShearMap();
  // Apply the mask on the shear map and perform inpainting on it
  double *gammaCorrArray = new double[Xaxis*Yaxis*Zaxis];
  for (size_t i=0; i<Xaxis; i++)
  {
    for (size_t j=0; j<Yaxis; j++)
    {
      unsigned int globalIndex = j*Xaxis +i;
      gammaCorrArray[globalIndex] = gammaMap.getBinValue(i, j, 0)*(1-(*m_maskValues)[i][j][0])
                                   +shearMap.getBinValue(i, j, 0)*(*m_maskValues)[i][j][0];
      globalIndex += Xaxis*Yaxis;
      gammaCorrArray[globalIndex] = gammaMap.getBinValue(i, j, 1)*(1-(*m_maskValues)[i][j][1])
                                   +shearMap.getBinValue(i, j, 1)*(*m_maskValues)[i][j][1];
    }
  }
  CoordinateBound m_CB = shearMap.getCoordinateBound();
  ShearMap gammaMapCorr(gammaCorrArray, Xaxis, Yaxis, Zaxis, m_CB, shearMap.getNumberOfGalaxies());

  delete [] gammaCorrArray;
  gammaCorrArray = nullptr;

  ConvergenceMap kappaMapBack = gammaMapCorr.getConvMap();

  return kappaMapBack;
}

}  // namespace LE3_2D_MASS_WL_CARTESIAN
