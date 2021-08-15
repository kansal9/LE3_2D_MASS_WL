/**
 * @file src/lib/SphericalIO.cpp
 * @date 12/24/20
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

#include "LE3_2D_MASS_WL_SPHERICAL/SphericalIO.h"

using namespace Euclid::WeakLensing::TwoDMass;

static Elements::Logging logger = Elements::Logging::getLogger("SphericalIO");

namespace Euclid {
 namespace WeakLensing {
  namespace TwoDMass {
   namespace Spherical {

SphericalIO::SphericalIO() {}

SphericalIO::SphericalIO(LE3_2D_MASS_WL_SPHERICAL::SphericalParam &SphParam): m_SphParam(SphParam) {}

void SphericalIO::writePrimaryHeader(const RecordHdu &hdu) {

  std::vector<Record<boost::any>> records = {
    { "RADESYS", "", "", "Equatorial coordinate system" },
    { "EQUINOX", "", "", "Equinox of celestial coordinate system (e.g. 2000)" },
    { "DATE-OBS", "", "", "Start of observation" }, //format ‘yyyy-mm-ddThh:mm:ss.sss’
    { "DATE-END", "", "", "End of observation" },
    { "SHE_SOFT", "LENSMC", "", "Origin of Shear Catalog" },
    { "SHE_SEL", "FITCLASS=0", "", "Origin of Shear Catalog" },
    { "SOFTNAME", "LE3_2D_MASS_WL_KS", "", "Software used to create the product" },
    { "SOFTVERS", SWVersion, "", "Software version" },
    { "NITREDSH", m_SphParam.getNItReducedShear(), "", "Number of iterations for reduced shear" },
    { "STDREDSH", m_SphParam.getRSSigmaGauss(),"", "gaussian smoothing sigma for reduced shear" },
    { "NITINP", m_SphParam.getNInpaint(),"", "Number of iterations for inpainting" },
    { "NSCINP", m_SphParam.getNInpScales(), "Number of scales for inpainting"},
    { "VARPERSC", m_SphParam.getEqualVarPerScale(), "True if equal variance per scale forced in"},
    { "FBMODE", m_SphParam.getBmodesZeros(), "" , "True if B-mode forced to zero in the gaps"},
    { "DENTYPE", "GAUSSIAN","", "denoising type"},
    { "GAUSSSTD", m_SphParam.getSigmaGauss(), "", "Standard deviation of gaussian smoothing"},
    { "FDRVAL", m_SphParam.getThresholdFDR(),"", "false discovery rate threshold"},
    { "NZBINS", m_SphParam.getNbins(),"", "number of redshift bins"},
    { "BAL_BINS", m_SphParam.getBalancedBins(), "", "True if balanced bins are applied"},
    { "NRESAMPL", m_SphParam.getNResamples(), "", "number of resampling of input dictionary"},
    { "MCREAL", "", "", ""},
    { "MCSEED", "", "", ""},
  };
  hdu.writeRecords(records);
}

void SphericalIO::writeHdu(const RecordHdu &hdu, int Nside) {

  std::vector<Record<boost::any>> records = {
    { "PIXTYPE", "HEALPIX", "", "HEALPIX Pixelisation" },
    { "ORDERING", "RING", "", "Pixel ordering scheme, either RING or NESTED" },
    { "NSIDE", Nside, "", "Resolution parameter of HEALPIX" },
    { "FIRSTPIX", "0", "", "First pixel # (0 based)" },
    { "LASTPIX", "", "", "Last pixel # (0 based)" },
    { "INDXSCHM","", "", "Indexing: IMPLICIT or EXPLICIT" },
    { "OBJECT", "0", "", "Sky coverage, either FULLSKY or PARTIAL" },
    { "ZMIN", "", "", "" },
    { "ZMAX","", "", "" },
  };
  hdu.writeRecords(records);
}

/*bool SphericalIO::write_Map (const std::string& filename, Healpix_Map<double>& map){
 fitshandle handleC = fitshandle();
 // permissions to be checked (R_OK, W_OK, X_OK) or the existence test (F_OK)
 /*if (access( filename.c_str(), F_OK ) != -1 ) {
  //handleC.open(filename);
   fitsfile *ptr;
   int status = 0;
   fits_open_file(&ptr, filename.c_str(), READWRITE, &status);
 } else {
  handleC.create(filename);
 }*/
 /* handleC.create(filename);
  write_Healpix_map_to_fits(handleC, map, PLANCK_FLOAT64);  //PLANCK_INT64 (float), PLANCK_FLOAT64 (Double)
  handleC.set_key("RADECSYS", std::string("FK5"), "World Coordinate System ");
  handleC.set_key("TCTYP2", std::string("HPX"), "X coordinate type ");
  handleC.set_key("TCTYP3", std::string("HPY"), "Y coordinate type ");
  handleC.set_key("SOFTNAME", std::string("LE3_2D_MASS_WL_KS"), "Software used to create the product");
  handleC.set_key("PROJ", std::string("NONE"), "processing is performed on the sphere");

  handleC.close();
return true;
}*/

bool SphericalIO::write_Map (const std::string& filename, Healpix_Map<double>& map, const std::string &colname){
  arr<double> myarr = map.Map();
  std::vector<double> mapValues;
  myarr.copyTo(mapValues);
  if (colname.compare("GALCOUNT") == 0 ) {
    std::string GalExt = "GALCOUNT_SPHERE";
    m_SphParam.setExtName(GalExt);
  }
  //mapValues.erase( std::remove(mapValues.begin(), mapValues.end(), 0), mapValues.end()); //removing Zero values
  mapValues.shrink_to_fit();
  if (access( filename.c_str(), F_OK ) != -1 ) {
    MefFile SphFile(filename, MefFile::Permission::Edit);
    const auto col = generateColumn<double>(colname, mapValues);
    long nbHdu = SphFile.hduCount();
    int nh = nbHdu-1;
    const auto &ext = SphFile.access<BintableHdu>(nh);
    /*if (ext.hasColumn(colname)) {
      const auto &ext_new = SphFile.assignBintableExt(m_SphParam.getExtName(), col);
      writeHdu(ext_new);
    } else {*/
      ext.appendColumn(col);
   // }
  } else {
    MefFile SphFile(filename, MefFile::Permission::Overwrite);
    const auto &primary = SphFile.accessPrimary<>();
    writePrimaryHeader(primary);
    const auto col = generateColumn<double>(colname, mapValues);
    const auto &ext = SphFile.assignBintableExt(m_SphParam.getExtName(), col);
    writeHdu(ext, map.Nside());
  }

 return true;
}
   } /* namespace Spherical */
  } /* namespace TwoDMass */
 } /* namespace WeakLensing */
} /* namespace Euclid */
