/**
 * @file src/lib/ReconstructMR.cpp
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

#include "LE3_2D_MASS_WL_CARTESIAN/ReconstructMR.h"

using namespace Euclid::WeakLensing::TwoDMass;

static Elements::Logging logger = Elements::Logging::getLogger("ReconstructMR");

namespace LE3_2D_MASS_WL_CARTESIAN {

ReconstructMR::ReconstructMR() {
  double Val=0.;
  int i=0, j=0;

  C1 = sqrt(2./M_PI);
  C2 = sqrt(2.);
  Step = 0.01;
  Np = (int) (5. / Step + 1.5);
  TabHsGauss = new double [Np];
  TabHnGauss = new double [Np];

  std::fill_n(TabHsGauss, Np, 0.);
  std::fill_n(TabHnGauss, Np, 0.);

  while (Val < 5.) {
      TabHsGauss[i] =  grad_hs_sig1(Val);
      TabHnGauss[i++] = grad_hn_sig1(Val);
      Val = Val + Step;
  }
}

ReconstructMR::ReconstructMR(const ReconstructMR& copy): C1(copy.C1), C2(copy.C2), Np(copy.Np),
   TabHsGauss(copy.TabHsGauss), TabHnGauss(copy.TabHnGauss), Step(copy.Step){}

ReconstructMR::~ReconstructMR()  {
   if (TabHnGauss  != NULL) {
     delete [] TabHnGauss;
   }
   if ( TabHsGauss != NULL) {
     delete [] TabHsGauss;
   }
   Np = 0;
}

double ReconstructMR::grad_hn_sig1(double Val) {
    int Sign = (Val >= 0) ? 1: -1;
    double Coef = Val*Sign;
    Coef = Coef*std::erfc(Coef/C2) + C1 * (1-exp(-(Coef*Coef)/2.));
    return Coef*Sign;
}

double ReconstructMR::grad_hs_sig1(double Val) {
     int Sign = (Val >= 0) ? 1: -1;
     double Coef = Val*Sign;
     Coef = -Coef*std::erf(Coef/C2) + C1 * (1-exp(-(Coef*Coef)/2.));
     return Coef*Sign;
}

double ReconstructMR::grad_hs(double Val) {
    int Ind;
    double ValRet = 0.;
    if (Val >= 5.) {
     ValRet = TabHsGauss[Np-1] - Val;
    } else if (Val > 0) {
       Ind = (int) (Val / Step + 0.5);
       ValRet = TabHsGauss[Ind];
    }
    return ValRet;
}

double ReconstructMR::grad_hn(double Val) {
   int Ind;
   double ValRet = 0.;
   if (Val >= 5.) ValRet = TabHnGauss[Np-1];
   else if (Val > 0)
   {
       Ind = (int) (Val / Step + 0.5);
       ValRet = TabHnGauss[Ind];
   }
   return ValRet;
}

double ReconstructMR::filter (double CoefDat, double Alpha, double Sig) {
   double ValRet=0.;
   double Sigma=Sig;
   if (Sigma < FLOAT_EPSILON) {
     Sigma = FLOAT_EPSILON;
   }
   if (fabs(Alpha) < FLOAT_EPSILON) {
     ValRet = CoefDat;
   } else if (Alpha < 0) {
     ValRet = 0.;
   } else {
      double ValDat = fabs(CoefDat) / Sigma;
      double CoefMax = ValDat + 3;
      double Coef, CoefMin = 0.;
      int Iter=0;
      do {
        Iter++;
        double Diff = 0.;
        Coef = (CoefMin+CoefMax)/2.;
        Diff =  grad_hs(ValDat-Coef) + Alpha * grad_hn(Coef);
        if (Diff > 0.) {
          CoefMax = Coef;
        } else {
          CoefMin = Coef;
        }
       } while (( CoefMax-CoefMin > 0.001) && (Iter < 100));
       Coef = (CoefDat >= 0.) ? Coef : -Coef;
       ValRet = Coef*Sigma;
   }
   return ValRet;
 }

}  // namespace LE3_2D_MASS_WL_CARTESIAN
