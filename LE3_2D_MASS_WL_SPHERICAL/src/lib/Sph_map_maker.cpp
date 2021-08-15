/**
 * @file src/lib/Sph_map_maker.cpp
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

#include "LE3_2D_MASS_WL_SPHERICAL/Sph_map_maker.h"
#include <cmath>
//#define M_PI 3.141592653589793238462643383279502884197;

using namespace Euclid::WeakLensing::TwoDMass;

static Elements::Logging logger = Elements::Logging::getLogger("SphericalMapMaker");
namespace LE3_2D_MASS_WL_SPHERICAL {

Sph_map_maker::Sph_map_maker(LE3_2D_MASS_WL_SPHERICAL::SphericalParam &SphParam):npix(0), m_SphParam(SphParam) {}

//std::pair<Healpix_Map<double>, Healpix_Map<double> > Sph_map_maker::create_ShearMap
//                 (std::vector<std::vector<double> > &catData){
std::tuple<Healpix_Map<double>, Healpix_Map<double>, Healpix_Map<double> > Sph_map_maker::create_ShearMap
                                                                (std::vector<std::vector<double> > &catData) {
  int order = hb.nside2order(m_SphParam.getNside());
  Healpix_Map<double> g1_hmap(order, RING);
  Healpix_Map<double> g2_hmap(order, RING);
  g1_hmap.fill(0.);
  g2_hmap.fill(0.);
  //npix = 12 * pow(nside,2);
  npix = g1_hmap.Npix();
  logger.info()<<"npix: "<<npix;

  Healpix_Map<double> ngal_hp(order, RING);
  ngal_hp.fill(0.);

  int ngal=catData[0].size();
  logger.info()<<"ngal: "<<ngal;
  //logger.info()<<"catData[1][0]: "<<catData[1][0];

 for (int i= 0; i<ngal; i++) {
  double theta, phi;
 /* if (m_SphParam.getZMin() < catData[5][i] < m_SphParam.getZMax()){
   //coordinates (theta and phi) in radian
   theta = -M_PI/180.0 * catData[1][i] + M_PI*double(0.5);
   phi = M_PI/180.0*(catData[0][i]);
  } else {*/
   //coordinates (theta and phi) in radian
   double dec = catData[1][i];
   if (isnan(dec)) {
     dec = 0.;
   }
  // theta = -M_PI/180.0 * catData[1][i] + M_PI*double(0.5);
   theta = -M_PI/180.0 * dec + M_PI*double(0.5);
   //logger.info()<<"theta: "<< theta;
   double ra = catData[0][i];
   if (isnan(ra)) {
     ra = 0.;
   }
  // phi = M_PI/180.0*(catData[0][i]);
   phi = M_PI/180.0*(ra);
  //}

   pointing ptg = pointing(theta, phi);
   ptg.normalize();
   auto id_pix = g2_hmap.ang2pix(ptg);
   auto id = g1_hmap.ang2pix(ptg);
   g1_hmap[id] = g1_hmap[id] + catData[3][i];
   g2_hmap[id_pix] = g2_hmap[id_pix] + catData[4][i];
   ngal_hp[id_pix] = ngal_hp[id_pix] + 1.;
 }
 for (unsigned int id_pix=0; id_pix<npix; id_pix++){
   if (ngal_hp[id_pix] != 0){
     g1_hmap[id_pix] = g1_hmap[id_pix]/ngal_hp[id_pix];
     g2_hmap[id_pix] = g2_hmap[id_pix]/ngal_hp[id_pix];
   }
 }
 //free(ngal_hp);
 //return std::pair<Healpix_Map<double>, Healpix_Map<double> > (g1_hmap, g2_hmap);
 return std::make_tuple(g1_hmap, g2_hmap, ngal_hp);
}

} // namespace LE3_2D_MASS_WL_SPHERICAL
