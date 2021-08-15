/**
 * @file LE3_2D_MASS_WL_UTILITIES/ReadCatalog.h
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

#ifndef _LE3_2D_MASS_WL_UTILITIES_READCATALOG_H
#define _LE3_2D_MASS_WL_UTILITIES_READCATALOG_H

#include "ElementsKernel/Logging.h"
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <boost/filesystem.hpp>

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"
#include "ElementsKernel/ProgramHeaders.h"

namespace Euclid {
 namespace WeakLensing {
  namespace TwoDMass {

/**
 * @class ReadCatalog
 * @brief
 *
 */
class ReadCatalog {

public:

 /**
   Empty Constructor
  */
  ReadCatalog();

 /**
   @brief Constructor
  */
  ReadCatalog(std::string catType);

  /**
   * @brief Destructor
   */
  virtual ~ReadCatalog() = default;

 /**
  * @brief         method returns the shear catalog data
  * @param[in]     <workdir> It's the working directory
  * @param[in]     <catalogName> It's the filename of the catalog
  * @param[in]     <data> <std::vector<std::vector<double> > > data of catalog in vectors column-wise
 */
  void readShearCatalog(fs::path& workdir, fs::path& catalogName, std::vector < std::vector < double> >& data);

 /**
  * @brief    The parameter which returns the type of input shear catalog used
 */
  std::string getCatalogReadType();
private:
 /**
  * @brief    The parameter which stores the type of input shear catalog used
 */
  std::string m_catalogType;

};  // End of ReadCatalog class

} /* namespace TwoDMass */
} /* namespace WeakLensing */
} /* namespace Euclid */
#endif
