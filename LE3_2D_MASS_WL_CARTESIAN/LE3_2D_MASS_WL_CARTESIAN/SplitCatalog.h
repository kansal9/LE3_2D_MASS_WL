/**
 * @file LE3_2D_MASS_WL_CARTESIAN/SplitCatalog.h
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

#ifndef _LE3_2D_MASS_WL_CARTESIAN_SPLITCATALOG_H
#define _LE3_2D_MASS_WL_CARTESIAN_SPLITCATALOG_H

#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

#include <vector>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <numeric>

namespace LE3_2D_MASS_WL_CARTESIAN {

/**
 * @class SplitCatalog
 * @brief
 *
 */
class SplitCatalog {

public:

  /**
   * @brief     Constructor
   *  @param    Cat_Data <vector<vector<double> > >, Input catalog data
   *  @param    param <CartesianParam>, Parameters read from parameter file
   *  @param    catType <std::string>, Origin of input catalogue
   *  @return   set object of class
   */
  SplitCatalog(std::vector<std::vector<double> >& Cat_Data, const CartesianParam& param, std::string& catType);

  /**
   * @brief     Constructor
   *  @param    Cat_Data <vector<vector<double> > >, Input catalog data
   *  @param    param <SphericalParam>, Parameters read from parameter file
   *  @param    catType <std::string>, Origin of input catalogue
   *  @return   set object of class
   */
  SplitCatalog(std::vector<std::vector<double> >& Cat_Data,
                  const LE3_2D_MASS_WL_SPHERICAL::SphericalParam& Sphparam, std::string& catType);

  /**
   * @brief Destructor
   */
  virtual ~SplitCatalog() = default;

  /**
   *  @brief    Private method to split catalog into sub-catalogs based on the need
   *            and writes the sub_catalogs into FITS format file
   *  @param    <Filenames>, name of the subcatalogs
   *  @param    <datadir>, path to data directory to save sub-catalogues
   *  @return   true if catalogs written successfully
  */
 bool writeSubCatalogs(fs::path& datadir, std::vector<std::string>& Filenames);

  /**
   *  @brief    method to split catalog into sub-catalogs based on the number of bins
   *  @param    <datadir>, path to data directory to save sub-catalogues
   *  @param    <Filenames>, name of the subcatalogs
  */
 void getSplittedCatalogs(fs::path& datadir, std::vector<std::string>& Filenames);

  /**
   *  @brief    method to split catalog into equal binned sub-catalogs
   *  @param    <datadir>, path to data directory to save sub-catalogues
   *  @param    <Filenames>, name of the subcatalogs
  */
 void getEqualBinSplittedCatalogs(fs::path& datadir, std::vector<std::string>& Filenames);

  /**
   *  @brief    write sub catalogs
   *  @param    Filenames <vector<string> >, name of the subcatalog
   *  @param    Data <vector<double> >, data to write
   *  @param    colname <string>, column name of the respective data
   *  @return   true if column is written successfully in the respective catalog
  */
 bool saveSubCatalogs (const std::string& filename, std::vector<double>& Data, const std::string &colname);

  /**
   *  @brief    returns indices of selected galaxies according to z bin
   *  @param    indices <vector<int> >, indices of selected galaxies
   *  @param    t_zmin <double>, min z value
   *  @param    t_zmax <double>, max z value
  */
 void getIndices(std::vector < int >& Indices, double t_zmin, double t_zmax);

  /**
   *  @brief    returns column name as in Data model
   *  @param    columnNames <vector<string> >, column names to be written in subcatalog
  */
 //void getColname(std::vector < std::string >& columnNames);

private:

  /**
   *  @brief <m_cartesianParam>, CartesianParam object with catalog parameters
  */
 LE3_2D_MASS_WL_CARTESIAN::CartesianParam m_cartesianParam;
  /**
   *  @brief <m_SphericalParam>, SphericalParam object with catalog parameters
  */
 LE3_2D_MASS_WL_SPHERICAL::SphericalParam m_SphericalParam;
  /**
   *  @brief <InData>, Catalog Data
  */
 std::vector<std::vector<double> > InData;
  /**
   *  @brief <m_nbZBins>, Number of redshift bins in input parameter file
  */
 size_t m_nbZBins;

  /**
   *  @brief <zMin>, Min Redshift value based on input parameter file
   *         <zMax>, Max Redshift value based on input parameter file
  */
 double m_zMin, m_zMax;
  /**
   *  @brief <m_balancedBins>, balancedBins input parameter
  */
 long m_balancedBins;
  /**
   *  @brief <m_catType>, origin of input catalog (helps to re-write sub-catalogues)
  */
 std::string m_catType;

};  // End of SplitCatalog class

/**
 *  @brief    Template to split a vector
 *  @param    vec <vector<T> >, vector to split
 *  @param    n <uint64_t>, number of sub-vectors
 *  @return   splitted vector
*/

template<typename T>
std::vector<T> split(std::vector<T>& vec, uint64_t n, size_t i) {
  std::vector<T> out_vec;
  uint64_t quotient = vec.size() / n;
  uint64_t reminder = vec.size() % n;

  auto start_itr = std::next(vec.cbegin(), i*quotient);;
 // auto end_itr2 = (reminder > 0) ? (quotient + !!(reminder--)) : quotient;
  auto end_itr = std::next(vec.cbegin(), i*quotient + quotient);
 // auto start_itr2 = i*quotient;

  //out_vec.push_back(std::vector<T>(vec.begin() + start_itr2, vec.begin() + end_itr2));

  // Sol1
  // code to handle the last sub-vector as it might
  // contain less elements
  if ((i*quotient + quotient) > vec.size()) {
      end_itr = vec.cend();
  }
  std::copy(start_itr, end_itr, back_inserter(out_vec));

//Sol2
 /*   if (i < reminder) {
        auto end_itr = std::next(vec.cbegin(), i*quotient + quotient + 1);
       std::copy(start_itr, end_itr, back_inserter(out_vec));
   }
    else if (i != n - 1) {
     auto end_itr = std::next(vec.cbegin(), i*quotient + quotient);
     std::copy(start_itr, end_itr, back_inserter(out_vec));
  }
    else {
     auto end_itr = vec.cend();
     std::copy(start_itr, end_itr, back_inserter(out_vec));
    }*/

return out_vec;
}

/**
 *  @brief  Extract elements from vector at indices specified by range
*/

template <typename ForwardIt, typename IndicesForwardIt>
inline std::vector<typename std::iterator_traits<ForwardIt>::value_type>
extract_at(
    ForwardIt first,
    IndicesForwardIt indices_first,
    IndicesForwardIt indices_last)
{
    typedef std::vector<typename std::iterator_traits<ForwardIt>::value_type>
        vector_type;
    vector_type extracted;
    extracted.reserve(static_cast<typename vector_type::size_type>(
        std::distance(indices_first, indices_last)));
    for(; indices_first != indices_last; ++indices_first)
        extracted.push_back(*(first + *indices_first));
    return extracted;
}

/**
 *  @brief  Extract elements from collection specified by collection of indices
*/

template <typename TVector, typename TIndicesVector>
inline TVector extract_at(const TVector& data, const TIndicesVector& indices)
{
    return extract_at(data.begin(), indices.begin(), indices.end());
}

}  // namespace LE3_2D_MASS_WL_CARTESIAN

#endif
