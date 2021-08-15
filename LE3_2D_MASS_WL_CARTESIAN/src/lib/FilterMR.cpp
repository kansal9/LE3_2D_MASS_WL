/**
 * @file src/lib/FilterMR.cpp
 * @date 03/12/21
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

#include "LE3_2D_MASS_WL_CARTESIAN/FilterMR.h"

using namespace Euclid::WeakLensing::TwoDMass;
using LE3_2D_MASS_WL_CARTESIAN::GetMap;
using LE3_2D_MASS_WL_CARTESIAN::Matrix;
using LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap;
using LE3_2D_MASS_WL_CARTESIAN::ReconstructMR;

static Elements::Logging logger = Elements::Logging::getLogger("FilterMR");

namespace LE3_2D_MASS_WL_CARTESIAN {

//FilterMR::FilterMR(ConvergenceMap &convMap, LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam):
//         m_cartesianParam(cartesianParam), convMap(convMap), m_MP(convMap.getXdim(), convMap.getYdim()) {

FilterMR::FilterMR(ConvergenceMap &convMap, bool PositiveCons, bool KillLastScale, bool KillIsol,
              int niter, int FirstScale, LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam):
              m_cartesianParam(cartesianParam), convMap(convMap), m_MP(convMap.getXdim(), convMap.getYdim()),
              m_KillLastScale(KillLastScale), m_KillIsol(KillIsol) {

 // Retrieve the size of the maps
 m_Xaxis = convMap.getXdim();
 m_Yaxis = convMap.getYdim();
 m_Zaxis = convMap.getZdim();
 m_PositveConstraint = PositiveCons;
 m_FirstScale = FirstScale;
 m_niter = niter;
 MaxIter = 20;
 m_SigmaNoise = 0.0;
 avar = 2.0;
 logger.info()<<"axis dim: "<<m_Xaxis<<" "<<m_Yaxis<<" "<<m_Zaxis;
 logger.info()<< std::boolalpha << "Kill Last Scale = " << m_KillLastScale;
 logger.info()<< std::boolalpha << "kill isolated pixels = " << m_KillIsol;
 logger.info()<< std::boolalpha << "apply positive constraint = " << m_PositveConstraint;
 logger.info()<< "First Scale = " << m_FirstScale - 1;
 m_FDR = m_cartesianParam.getThreshold();
 if (m_FDR == 0) {
   m_FDR = 0.05;
 }
 logger.info()<< "FDR_val = " << m_FDR;
 nbScales = m_cartesianParam.getnbScales();
 if (nbScales==0) {
  nbScales = int(log(m_Xaxis)/log(2.))-2.;
 }
 logger.info()<<"nbScales: " << nbScales;
}

FilterMR::FilterMR(const FilterMR& copy): convMap(copy.convMap),
  m_MP(copy.m_MP), m_Xaxis(copy.m_Xaxis), m_Yaxis(copy.m_Yaxis), m_Zaxis(copy.m_Zaxis),
  m_FirstScale(copy.m_FirstScale), m_SigmaNoise(copy.m_SigmaNoise), nbScales(copy.nbScales), m_FDR(copy.m_FDR),
  m_PositveConstraint(copy.m_PositveConstraint), m_KillLastScale(copy.m_KillLastScale), m_KillIsol(copy.m_KillIsol),
  m_niter(copy.m_niter), avar(copy.avar) {}


double FilterMR::fdr_pv(std::vector<double> &PValue, double alpha) {
   std::vector<double> input_pval;
    for (size_t i=0; i < PValue.size(); i++){
        input_pval.push_back( PValue[i]);
    }
    std::sort(input_pval.begin(), input_pval.end());
//    input_pval[1] =0;
    int m = input_pval.size();
    double p_cutoff=0.;
    for (int i = 0; i < m; i++) {
      double a = alpha * i / m;
      if (input_pval[i] < a) {
        p_cutoff = input_pval[i];
      }
   }
   return p_cutoff;
}

double FilterMR::myErfInv(double x){
   double tt1, tt2, lnx, sgn;
   sgn = (x < 0) ? -1.0 : 1.0;
   x = (1 - x)*(1 + x); // x = 1 - x*x;
   lnx = log(x);

   tt1 = 2/(M_PI*0.147) + 0.5 * lnx;
   tt2 = 1/(0.147) * lnx;

   return (sgn*sqrt(-tt1 + sqrt(tt1*tt1 - tt2)));
}

double FilterMR::xerfc (double F) {
   double Nu;
   double P = F;
   if (P >= 1.) {
     P = 1.;
   } else {
    if (P < 0) {
       P = 0.;
    }
    if (P > 0.5) {
       Nu = sqrt(2.)*myErfInv( 1 - ((1. - P)/0.5));
    } else {
       Nu = - sqrt(2.)*xerfc((double) P/0.5);
    }
   }
   return Nu;
}

double FilterMR::getProbility(double val, double sigLevel) {
 double P=0., vn=0.;
 if (fabs(val)<FLOAT_EPSILON) {
   P = 1.;
 } else {
  if (sigLevel<FLOAT_EPSILON) {
    P = 0.;
  }
  else {
   vn = fabs(val)/(sqrt(2.)*sigLevel);
   if (vn>3.5) {
     P = 0.;
   } else {
     P = std::erfc(vn);
   }
  }
 }
//std::cout << "P:" << P << std::endl;
 return P;
}

void FilterMR::get_fdr(std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix>& myBand, std::vector<double>& NSigma) {
 logger.info() << "Begin";
 
 for (int kScale = 0; kScale<nbScales-1; kScale++){
   std::vector < double> PValue;
  // logger.info()<< "Sig " << m_SigmaNoise * Euclid::WeakLensing::TwoDMass::NormB3Spline[kScale];
    for (int j=0; j<m_Yaxis; j++){
       for (int i=0; i<m_Xaxis; i++){
         PValue.push_back(getProbility(myBand[kScale].getValue(i, j),
                (m_SigmaNoise * Euclid::WeakLensing::TwoDMass::NormB3Spline[kScale])));
       }
     }
    double alpha_b = m_FDR * pow(avar, double(kScale));
    //alpha_b = (kScale==0) ? 1. - erf(NSigma[kScale] / sqrt(2.)): alpha_b*2;
    if (alpha_b > 0.5) {
     alpha_b = 0.5;
    }
    double Pdet = fdr_pv(PValue, alpha_b);

    double temp = 0.5+(1. - Pdet)/2.;
   NSigma[kScale] = fabs(xerfc(temp));
   if ((NSigma[kScale] < 7) && (NSigma[kScale] > 0)) {
     NSigma[kScale] = fabs(xerfc(temp));
   } else {
     NSigma[kScale] = 7.;
   }

   logger.info()<<"band: " << kScale + 1 <<"    m_FDR: " << alpha_b <<"    nsigma: " << NSigma[kScale];
 }

 logger.info() << "End";
}

void FilterMR::kill_isolated(LE3_2D_MASS_WL_CARTESIAN::Matrix& image) {
  LE3_2D_MASS_WL_CARTESIAN::Matrix has_neighbor(m_Xaxis, m_Yaxis);

   has_neighbor.reset(); // set all values to 0
 /*  for (int j=1; j<m_Yaxis-1; j++){
     for (int i=1; i<m_Xaxis-1; i++){
        if (image.getValue(i,j)-image.getValue(i-1,j) > 0 && image.getValue(i,j)-image.getValue(i+1,j) > 0
            && image.getValue(i,j)-image.getValue(i,j-1) > 0 && image.getValue(i,j)-image.getValue(i,j+1) > 0) {
           has_neighbor.setValue(i, j, 1);
        }      
     }
   }*/
   for (int j=1; j<m_Yaxis-1; j++){
     for (int i=1; i<m_Xaxis-1; i++){
        if (has_neighbor.getValue(i,j) || image.getValue(i-1,j) > 0 && has_neighbor.getValue(i,j) ||
            image.getValue(i+1,j) > 0 && has_neighbor.getValue(i,j) || image.getValue(i,j-1) > 0 &&
            has_neighbor.getValue(i,j) || image.getValue(i,j+1) > 0) {
           has_neighbor.setValue(i, j, 1);
        }      
     }
   }
   for (int j=0; j<m_Yaxis; j++){
     for (int i=0; i<m_Xaxis; i++){
      if (has_neighbor.getValue(i,j) == 1) {
         image.setValue(i, j, 0);
      }
     }
   }
}

void FilterMR::set_support (std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix> & myBand,
                            std::vector<double>& NSigma, int it) {

    for (int j=0; j<m_Yaxis; j++){
       for (int i=0; i<m_Xaxis; i++){
         double temp = myBand[it].getValue(i, j);
         if (temp < (NSigma[it]*m_SigmaNoise*Euclid::WeakLensing::TwoDMass::NormB3Spline[it])) {
           myBand[it].setValue(i, j, 0.);
         }
         if (temp > 0) {
           myBand[it].setValue(i, j, 1.);
         }
       }
     }
}

void FilterMR::mr_support (std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix>& myBand, std::vector<double>& NSigma,
                      std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix>& mr_myBand) {

  std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix> copy_myBand;
  copy_myBand.assign(myBand.begin(), myBand.end());
  for (int it = m_FirstScale-1; it < nbScales-1; it++) {
     set_support(copy_myBand, NSigma, it);

     if (m_KillIsol == true) {
       kill_isolated(copy_myBand[it]);
     }
     for (int j=0; j<m_Yaxis; j++){
      for (int i=0; i<m_Xaxis; i++){
        mr_myBand[it].setValue(i, j, copy_myBand[it].getValue(i, j));
      }
     }
  }
}

ConvergenceMap* FilterMR::performFiltering() {
 // Create a copy of the conv map
 ConvergenceMap *kappaMap = new ConvergenceMap(convMap);

 // Allocate memory for kappa
  Matrix kappaE (m_Xaxis, m_Yaxis);
  Matrix kappaB (m_Xaxis, m_Yaxis);

 for (int j=0; j<m_Yaxis; j++){
  for (int i=0; i<m_Xaxis; i++){
   kappaE.setValue(i, j, kappaMap->getBinValue(i, j, 0));
   kappaB.setValue(i, j, kappaMap->getBinValue(i, j, 1));
  }
 }

 m_SigmaNoise = kappaE.getStandardDeviation();
 logger.info()<<"sigmaNoise: " << m_SigmaNoise;

 std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix> myBand = m_MP.transformBspline(kappaE, nbScales);

 std::vector<double> NSigma_l(nbScales-1);

 get_fdr(myBand, NSigma_l);

 std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix> mr_myBand = myBand;
 for (size_t b = 0; b < mr_myBand.size(); b++) {
     mr_myBand[b].reset();
 }

 mr_support(myBand, NSigma_l, mr_myBand);

 mr_myBand[nbScales-1].reset();

 applyfilter(myBand, NSigma_l);

 kappaE.reset();

 kappaE = reconstruct_optimized (myBand, mr_myBand);

 if (true == m_PositveConstraint) {
  kappaE.applyThreshold(0.);
 }

 delete kappaMap;

 // Allocate memory for the kappa array
 double *kappaArray = new double[m_Xaxis*m_Yaxis*m_Zaxis];
 // Fill the convergence map array
 for (int i=0; i<m_Xaxis; i++){
  for (int j=0; j<m_Yaxis; j++) {
   kappaArray[j*m_Xaxis +i] = kappaE.getValue(i, j);
   //kappaArray[m_Xaxis*m_Yaxis + j*m_Xaxis +i] = 0;
  }
 }

 ConvergenceMap kappaMap1(kappaArray, m_Xaxis, m_Yaxis, m_Zaxis);

 kappaMap = new ConvergenceMap(kappaMap1);
 return kappaMap;
}

Matrix FilterMR::reconstruct_optimized(std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix> & band,
                  std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix> & mr_band) {

   for (int b = 0; b < m_FirstScale-1; b++) {
       band[b].reset();
   }
 
   //Kill last scale
   if (true == m_KillLastScale) {
          band[nbScales-1].reset();
   }

  Matrix map (m_Xaxis, m_Yaxis);
  map = reconstruct (band);

  std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix> band_i;
  std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix> band_res;

  for (int it =0; it < m_niter; it++) {
    std::cout << "start iteration " << it+1 << std::endl;
    band_i = m_MP.transformBspline(map, nbScales);
    band_res = band_i;
    for (int b = 0; b < nbScales; b++) {
     for (int j=0; j<m_Yaxis; j++){
      for (int i=0; i<m_Xaxis; i++){
          double a0 = band[b].getValue(i, j);
          double b0 = band_i[b].getValue(i, j);
          band_res[b].setValue(i, j, (a0-b0));
          double c0 = mr_band[b].getValue(i, j);
          if (c0 == 0) {
           band_res[b].setValue(i, j, 0.);
          }
      }
     }
    }

    for (int j=0; j<m_Yaxis; j++){
     for (int i=0; i<m_Xaxis; i++){
          double a1 = band[nbScales-1].getValue(i, j);
          double b1 = band_i[nbScales-1].getValue(i, j);
          band_res[nbScales-1].setValue(i, j, a1-b1);
     }
    }
    Matrix res(m_Xaxis, m_Yaxis);
    res = reconstruct (band_res);

    for (int j=0; j<m_Yaxis; j++){
     for (int i=0; i<m_Xaxis; i++){
      double a2 = map.getValue(i, j);
      double b2 = res.getValue(i, j);
      if (b2>0) {
        map.setValue(i, j, a2+b2);
      }
     }
    }
  }
 return map;
}

Matrix FilterMR::reconstruct(std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix> & band) {
 Matrix Image (m_Xaxis, m_Yaxis);
 Matrix temp (m_Xaxis, m_Yaxis);
 for (int j=0; j<m_Yaxis; j++){
   for (int i=0; i<m_Xaxis; i++){
     Image.setValue(i, j, band[nbScales-1].getValue(i, j));
   }
 }

 size_t s = nbScales-2;
 do {
   temp = m_MP.smoothBspline(Image, s);
  for (int j=0; j<m_Yaxis; j++){
   for (int i=0; i<m_Xaxis; i++){
        double aa = temp.getValue(i, j);
        double bb = band[s].getValue(i, j);
        Image.setValue(i, j, aa+bb);
     }
   }

 } while ( s-- );

 return Image;
}

void FilterMR::applyfilter(std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix>& band, std::vector<double>& NSigma) {

 std::vector<double> TabAlpha (nbScales-1);
 std::fill (TabAlpha.begin(), TabAlpha.end(), 1.);

 TabAlpha[0] = 5.;

 std::vector<double> TabDelta(nbScales-1);

 std::vector<double> regulMin(nbScales-1);
 std::fill (regulMin.begin(), regulMin.end(), 0.);
 std::vector<double> regulMax(nbScales-1);
 std::fill (regulMax.begin(), regulMax.end(), 100.);

 int iter=0;

 std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix> copyBand;// = band;
  //copy(band.begin(), band.end()-1, back_inserter(copyBand))

 ReconstructMR restore;

 do {
  iter++;
  logger.info() << " START Multiscale Entropy Iteration = " << iter;
  copyBand.assign(band.begin(), band.end());

  for(int b=0; b<nbScales-1; b++) {
    TabAlpha[b] = (regulMin[b] + regulMax[b])/2.;
    //logger.info() << " TabAlpha = " << TabAlpha[b];
  }

  fixedAlphaFilter(copyBand, TabAlpha, NSigma, restore);

  //copyBand[nbScales-1].reset();

  LE3_2D_MASS_WL_CARTESIAN::Matrix image (m_Xaxis, m_Yaxis);
  image = m_MP.reconsBspline(copyBand);

  double Flux=image.getFlux();
  double Mean=Flux / (m_Yaxis*m_Xaxis);

  logger.info() << " Min = " << image.getMin() << "    Max = " << image.getMax() << "    Flux = " << Flux;
  logger.info() << "    Mean = " << Mean << "     Sigma = " << image.getStandardDeviation();

  if (true == m_PositveConstraint) {
   image.applyThreshold(0.);
  }

   copyBand = m_MP.transformBspline(image, nbScales);
   estimateNewAlphaParam(band, copyBand, regulMin, regulMax, TabAlpha, TabDelta);
  logger.info() << " END Multiscale Entropy Iteration = " << iter;
  //std::cout << "max val: " << (fabs(*max_element(TabDelta.begin(), TabDelta.end())))<< std::endl;
 } while(((fabs(*max_element(TabDelta.begin(), TabDelta.end()))) > 0.01) && (iter < MaxIter));

  for (int b = 0; b<nbScales-1; b++){
   TabAlpha[b] = (regulMin[b] + regulMax[b])/2;
   //logger.info() << " TabAlpha = " << TabAlpha[b];
  }

  fixedAlphaFilter(band, TabAlpha, NSigma, restore);

}

void FilterMR::estimateNewAlphaParam(std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix> & band,
     std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix> & band_sol, std::vector<double>& regulMin,
     std::vector<double>& regulMax, std::vector<double>& TabAlpha, std::vector<double>& TabDelta) {
   double sigmaNoise;
  // for (int b=0; b<nbScales-1; b++) {
   for (int b=m_FirstScale-1; b<nbScales-1; b++) {
     sigmaNoise = 0.;
     double sigmaCoeff = m_SigmaNoise * Euclid::WeakLensing::TwoDMass::NormB3Spline[b];
     //double sigmaCoeff = Euclid::WeakLensing::TwoDMass::NormB3Spline[b];
     for (int j=0; j<m_Yaxis; j++){
      for (int i=0; i<m_Xaxis; i++){
       double residual = band[b].getValue(j, i) - band_sol[b].getValue(j, i);
       double var = sigmaCoeff*sigmaCoeff;
 	   sigmaNoise = sigmaNoise + ((residual * residual) / var);
       //if (var > FLOAT_EPSILON) {
 	      //sigmaNoise = sigmaNoise + ((residual * residual) / var);
       //} else {
       //   sigmaNoise = sigmaNoise + 1. + (residual * residual);
      // }
      }
     }
     sigmaNoise = sqrt(sigmaNoise/(double) (m_Yaxis*m_Xaxis));
	 if (sigmaNoise >= 1) {
       regulMax[b] =  TabAlpha[b];
     } else {
       regulMin[b] = TabAlpha[b];
     }
	 TabDelta[b] = fabs(sigmaNoise - 1.);
     //logger.info() << "scale = " << b << "     Normalized sig = " << sigmaNoise;
     //logger.info() << "Alpha = " << TabAlpha[b] << "     Delta = " << TabDelta[b];
   }
}

void FilterMR::fixedAlphaFilter(std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix>& band, std::vector<double>& TabAlpha,
                 std::vector<double>& NSigma, LE3_2D_MASS_WL_CARTESIAN::ReconstructMR & restore) {

 double alphaP, val=0.;
 for (int b = 0; b < nbScales-1; b++) {
   //  std::cout << "rp: " << TabAlpha[b] <<std::endl;
   double sigma = m_SigmaNoise * Euclid::WeakLensing::TwoDMass::NormB3Spline[b];
   for (int j=0; j<m_Yaxis; j++){
     for (int i=0; i<m_Xaxis; i++){
         if (b < m_FirstScale-1) {
           band[b].setValue(i, j, 0.);
         } else {
           double rp = TabAlpha[b];
           val = band[b].getValue(i, j);
           if (val>=0) {
             alphaP = fabs(val / (sigma * NSigma[b]));
             if (alphaP > 1) {
               alphaP = 1;
             }
             if (alphaP < FLOAT_EPSILON) {
               rp = -1;
             }
	         if (alphaP != 0) {
               rp *= (1-alphaP)/alphaP;
             }
           }
          // std::cout << "regularParam: " << rp << std::endl;
           double temp_val = restore.filter(val, rp, sigma);
           if (rp > 0) {
             band[b].setValue(i, j, temp_val);
           } else {
             band[b].setValue(i, j, val);
           }
         }
     }
   }
 }
}

}  // namespace LE3_2D_MASS_WL_CARTESIAN
