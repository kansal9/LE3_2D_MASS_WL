/**
 * @file src/program/LE3_2D_MASS_WL_PatchesToSphere.cpp
 * @date 02/21/20
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

#include <map>
#include <string>

#include <boost/program_options.hpp>
#include "ElementsKernel/ProgramHeaders.h"

#include "LE3_2D_MASS_WL_CARTESIAN/GetMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ShearMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/Projection.h"
#include <algorithm>

#include <cmath>
#include <vector>
#include "LE3_2D_MASS_WL_CARTESIAN/MatrixProcess.h"

#include "fitsio.h"

#include "healpix_cxx/fitshandle.h"
#include "healpix_cxx/pointing.h"
#include "healpix_cxx/healpix_map_fitsio.h"
#include "healpix_cxx/healpix_base.h"
#include "healpix_cxx/healpix_map.h"
#include "healpix_cxx/healpix_data_io.h"
#include "healpix_cxx/healpix_map_fitsio.h"

#include <healpix_cxx/alm.h>
#include <healpix_cxx/alm_fitsio.h>
#include <healpix_cxx/alm_healpix_tools.h>

using boost::program_options::options_description;
using boost::program_options::variable_value;

using std::string;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

using namespace Euclid::WeakLensing::TwoDMass;
using LE3_2D_MASS_WL_CARTESIAN::ShearMap;
using LE3_2D_MASS_WL_CARTESIAN::Projection;
using LE3_2D_MASS_WL_CARTESIAN::MatrixProcess;

static Elements::Logging logger = Elements::Logging::getLogger("LE3_2D_MASS_WL_PatchesToSphere");

class LE3_2D_MASS_WL_PatchesToSphere : public Elements::Program {

public:

  options_description defineSpecificProgramOptions() override {
  
    options_description options {};
   // input working directory
   options.add_options()
   ("workdir", po::value<string>()->default_value(""), "Work Directory");

   // input log directory: default is none
   options.add_options()
   ("logdir", po::value<string>()->default_value(""), "logs Directory");

   // input catalog file
   options.add_options()
   ("input", po::value<string>()->default_value(""), "input maps");

    //
    // !!! Implement the program options here !!!
    //
    return options;
  }

  Elements::ExitCode mainMethod(std::map<std::string, variable_value>& args) override {

    Elements::Logging logger = Elements::Logging::getLogger("LE3_2D_MASS_WL_PatchesToSphere");

    logger.info("#");
    logger.info("# Entering mainMethod()");
    logger.info("#");

    int nside = 512;
////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Get the workdir
////////////////////////////////////////////////////////////////////////////////////////////////////////
   fs::path workdir {args["workdir"].as<string>()};
   fs::path logPath (args["logdir"].as<string>());
   fs::path input {args["input"].as<std::string>()};

    ShearMap *m_shear;
    m_shear = new ShearMap((workdir/input).native());
  int x = m_shear->getXdim();
  int y = m_shear->getYdim();
  int z = m_shear->getZdim();
  int raMin = 45;
  int raMax = 55;
  int decMin = 15;
  int decMax = 25;

 //create object to get correct projection
 Projection projection;
 // Get the map characteristics
 double ra0 = 0.5*(raMin + raMax); // in degree
 double dec0 = 0.5*(decMin + decMax); // in degree
 double raRange = (raMax - raMin) * M_PI/180.; //in rad
 double decRange = (decMax - decMin) * M_PI/180.; //in rad    
  LE3_2D_MASS_WL_CARTESIAN::Matrix map(x, y);
  for (unsigned int j=0; j<x; j++) {
    for (unsigned int i=0; i<y; i++) {
      map.setValue(i, j, m_shear->getBinValue(i, j, 0));
    }
  }
 std::vector<double> ra;
 std::vector<double> dec;
 std::vector<double> data;
  // Loop over all pixels
  for (unsigned int i=0; i<x; i++) {
   for (unsigned int j=0; j<y; j++) {

       // Perform transform from pixel location to ra and dec
       double tmpx = (i+0.5)*raRange/x-0.5*raRange;
       double tmpy = (j+0.5)*decRange/y-0.5*decRange;
       std::pair<double, double> radec = projection.getInverseGnomonicProjection(tmpx, tmpy, ra0, dec0);
       ra.push_back(radec.first);
       dec.push_back(radec.second);
       data.push_back(map.getValue(i, j));
   }
 }	

Healpix_Base hb;
int order = hb.nside2order(nside);
 Healpix_Map<double> g_hmap(order, RING);
 g_hmap.fill(0.);
 int npix = g_hmap.Npix();
 logger.info()<<"npix: "<<npix;
 int *ngal_hp;
 //int *ngal_hp = new int[npix];
 ngal_hp = (int *)calloc(npix, sizeof(int));

 //initialize
 for (unsigned int n=0; n<npix; n++){
 ngal_hp[n]=0;
 }

 int ngal=ra.size();
 logger.info()<<"ngal: "<<ngal;

for (unsigned int i= 0; i<ngal; i++) {
 double theta, phi;
  theta = -M_PI/180.0 * dec[i] + M_PI*double(0.5);
  phi = M_PI/180.0*(ra[i]);

  pointing ptg = pointing(theta, phi);
  ptg.normalize();
  auto id_pix = g_hmap.ang2pix(ptg);
  g_hmap[id_pix] = g_hmap[id_pix] + data[i];
  ngal_hp[id_pix] = ngal_hp[id_pix] + 1.;
}
for (unsigned int id_pix=0; id_pix<npix; id_pix++){
  if (ngal_hp[id_pix] != 0){
    g_hmap[id_pix] = g_hmap[id_pix]/ngal_hp[id_pix];
  }
}
free(ngal_hp);

  fitshandle handleC = fitshandle();
  handleC.create("testSphere.fits");
  write_Healpix_map_to_fits(handleC, g_hmap, PLANCK_FLOAT64);  //PLANCK_INT64 (float), PLANCK_FLOAT64 (Double)
  handleC.set_key("RADECSYS", std::string("FK5"), "World Coordinate System ");
  handleC.set_key("TCTYP2", std::string("HPX"), "X coordinate type ");
  handleC.set_key("TCTYP3", std::string("HPY"), "Y coordinate type ");
  handleC.close();

    //
    // !!! Implement you program here !!!
    //

    logger.info("#");
    logger.info("# Exiting mainMethod()");
    logger.info("#");

    return Elements::ExitCode::OK;
  }

};

MAIN_FOR(LE3_2D_MASS_WL_PatchesToSphere)
