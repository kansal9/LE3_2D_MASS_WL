/**
 * @file src/lib/ReadCatalog.cpp
 * @date 05/15/21
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

#include "LE3_2D_MASS_WL_UTILITIES/ReadCatalog.h"

namespace fs = boost::filesystem;
using namespace Euclid::WeakLensing::TwoDMass;

static Elements::Logging logger = Elements::Logging::getLogger("ReadCatalog");

namespace Euclid {
 namespace WeakLensing {
  namespace TwoDMass {

ReadCatalog::ReadCatalog(): m_catalogType("LENSMC") {}
ReadCatalog::ReadCatalog(std::string catType): m_catalogType(catType) {}

void ReadCatalog::readShearCatalog(fs::path& workdir, fs::path& catalogName,
                                 std::vector < std::vector < double> >& data) {
   if (true == catalogName.string().empty()) {
     throw Elements::Exception() << "Input catalogue file name is not found . . . ";
   }
   fs::path datadir {workdir / "data"};
  // Variable to save input catalog name
   std::string inputCatalog;
   //std::vector<std::vector<double> > data;
   CatalogData galdata;
// check whether input shear catalog file is in fits format(standalone) / XML(pipeline) and then fetch Catalogfile name
   if (true == checkFileType(datadir /catalogName, Euclid::WeakLensing::TwoDMass::signFITS)) {
    if (false == fs::exists(datadir/catalogName)) {
       throw Elements::Exception() << "Input data product " << datadir/catalogName << " not found";
    }
     logger.info("Input Shear Catalog is in Fits format..");
     inputCatalog = (datadir /catalogName).string();
     galdata.readCatalog(inputCatalog, data);
     logger.info("Done reading Input Shear Catalog");
   } else { //Case 2: When input catalog file is in XML format
     if (false == fs::exists(workdir/catalogName)) {
       throw Elements::Exception() << "Input data product " << workdir/catalogName << " not found";
     }
     logger.info("Input Shear Catalog is in XML format..");
     galdata.getCatalogData(workdir, catalogName, data);
     logger.info("Done reading Input Shear Catalog");
   }
     m_catalogType = galdata.getCatalogType();
 //return data;
}

std::string ReadCatalog::getCatalogReadType() {
  return m_catalogType;
}

} /* namespace TwoDMass */
} /* namespace WeakLensing */
} /* namespace Euclid */
