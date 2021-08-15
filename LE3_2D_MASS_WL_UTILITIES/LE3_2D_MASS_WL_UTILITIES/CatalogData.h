/**
 * @file LE3_2D_MASS_WL_UTILITIES/CatalogData.h
 * @date 10/13/20
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

#ifndef _CATALOGDATA_H
#define _CATALOGDATA_H

#include <vector>
#include <cstdio>
#include <iostream>
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"
#include <boost/filesystem.hpp>
#include "ElementsKernel/Logging.h"

#include "EL_FitsFile/MefFile.h"

namespace fs = boost::filesystem;

namespace Euclid {
 namespace WeakLensing {
  namespace TwoDMass {

/**
 * @class CatalogData
 * @brief
 *
 */
class CatalogData {

public:

  /**
   * @brief Destructor
   */
  virtual ~CatalogData() = default;
 /**
   Empty Constructor
  */
  CatalogData();
 /**
   @brief Constructor
  */
  CatalogData(std::string catType);

 /**
  * @brief         check catalog format (fits/xml) and return data of that catalog
  * @param[in]     <workdir> <boost::filesystem::path> work directory
  * @param[in]     <InCatalog> <boost::filesystem::path> Input Catalog name
  * @param[in]     <data> <std::vector<std::vector<double> > > data of catalog in vectors column-wise
 */
  void getCatalogData(fs::path& workdir, fs::path& InCatalog, std::vector<std::vector<double> >& data);
  //std::vector<std::vector<double> > getCatalogData(fs::path& workdir, fs::path& InCatalog);

 /**
  * @brief         method to read the catalog and return data
  * @param[in]     <filename> It's the filename of the catalog with path
  * @param[in]     <Cat_data> <std::vector<std::vector<double> > > data of catalog in vectors column-wise
 */
  //std::vector<std::vector<double> > readCatalog(const std::string& filename);
  void readCatalog(const std::string& filename, std::vector<std::vector<double> >& Cat_Data);

 /**
  * @brief         method to read the cluster catalog and return data
  * @param[in]     <filename> It's the filename of the catalog
  * @param[in]     <Cat_data> <std::vector<std::vector<double> > > data of catalog in vectors column-wise
 */
  void readClusterCatalog(const std::string& filename, std::vector<std::vector<double> >& Cat_Data);

 /**
  * @brief         method to read the shear catalog and return data
  * @param[in]     <filename> It's the filename of the catalog
  * @param[in]     <Cat_data> <std::vector<std::vector<double> > > data of catalog in vectors column-wise
 */
  void readShearCatalog(const std::string& filename, std::vector<std::vector<double> >& Cat_Data);

 /**
  * @brief    The parameter which returns the type of input shear catalog used
 */
  std::string getCatalogType();

 /**
  * @brief    The parameter which returns the name of fits input shear catalog
 */
  std::string getCatalogFitsName();

private:
 /**
  * @brief    The parameter which stores the type of input shear catalog used
 */
  std::string m_catType;
 /**
  * @brief    The parameter which stores the input catalog name
 */
  std::string m_inputCatalog;

};  // End of CatalogData class
} /* namespace TwoDMass */
} /* namespace WeakLensing */
} /* namespace Euclid */
#endif
