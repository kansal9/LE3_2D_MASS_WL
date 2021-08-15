/**
 * @file LE3_2D_MASS_WL_SPHERICAL/SphericalIO.h
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

#ifndef _SPHERICALIO_H
#define _SPHERICALIO_H

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "ElementsKernel/Logging.h"

#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"

namespace Euclid {
 namespace WeakLensing {
  namespace TwoDMass {
   namespace Spherical {

/**
 * @class SphericalIO
 * @brief
 *
 */
class SphericalIO {

public:

  /**
   * @brief empty Constructor
   */
  SphericalIO();

  /**
   * @brief Constructor
   */
  SphericalIO(LE3_2D_MASS_WL_SPHERICAL::SphericalParam &SphParam);

  /**
   * @brief Destructor
   */
  virtual ~SphericalIO() = default;

  /**
  @brief  This method write Fits BinTable (according to DataModel) using input healpix map and filename
  @param  <filename> name of the map in fits format
  @param  <map> healpix map
  @param  <colname> column name in table
  @retun  true if map is wriiten correctly
  */
  bool write_Map (const std::string& filename, Healpix_Map<double>& map, const std::string &colname);

  /**
  @brief  This method write Map using input healpix map and filename
  @param  <filename> name of the map in fits foramt
  @param  <map> healpix map
  @retun  true if map is wriiten correctly
  */
  //bool write_Map (const std::string& filename, Healpix_Map<double>& map);

 /**
 * @brief   writes records to the Primary header
 * @param[in] Primary hdu to which records need to be written
 * This method writes records to the Primary header
 **/
  void writePrimaryHeader(const RecordHdu &hdu);

 /**
 * @brief   writes records to the given header
 * @param[in] hdu (other than primary) to which records need to be written
 * This method writes records to the given header
 **/
  void writeHdu(const RecordHdu &hdu, int Nside);

private:
  /**
   *  @brief <m_SphParam>, SphericalParam object with input parameters
  */
 LE3_2D_MASS_WL_SPHERICAL::SphericalParam m_SphParam;

};  // End of SphericalIO class

} /* namespace Spherical */
} /* namespace TwoDMass */
} /* namespace WeakLensing */
} /* namespace Euclid */
#endif
