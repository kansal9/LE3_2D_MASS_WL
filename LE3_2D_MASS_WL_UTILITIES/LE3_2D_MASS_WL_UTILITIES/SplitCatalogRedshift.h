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
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"

#include <vector>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <numeric>

namespace LE3_2D_MASS_WL_UTILITIES
{

/**
 * @class SplitCatalog
 * @brief
 *
 */
class SplitCatalogRedshift
{

public:

    /**
     * @brief    Constructor
     * @param    cat input CatalogData
     * @param    parameters
     * @param    catType <std::string>, origin of input catalogue
     */
    SplitCatalogRedshift(CatalogData& cat, const GenericParam& param,
                         std::string& catType);

    /**
     * @brief Destructor
     */
    virtual ~SplitCatalogRedshift() = default;

    /**
     *  @brief    Private method to split catalog into sub-catalogs based on the need
     *            and writes the sub_catalogs into FITS format file
     *  @param    <Filenames>, name of the subcatalogs
     *  @param    <datadir>, path to data directory to save sub-catalogues
     */
    void writeSubCatalogs(fs::path& datadir, std::vector<fs::path>& Filenames);

    /**
     *  @brief    method to split catalog into sub-catalogs based on the number of bins
     *  @param    <datadir>, path to data directory to save sub-catalogues
     *  @param    <Filenames>, name of the subcatalogs
     */
    void getSplittedCatalogs(fs::path& datadir,
                             std::vector<fs::path>& Filenames);

    /**
     *  @brief    method to split catalog into equal binned sub-catalogs
     *  @param    <datadir>, path to data directory to save sub-catalogues
     *  @param    <Filenames>, name of the subcatalogs
     */
    void getEqualBinSplittedCatalogs(fs::path& datadir,
                                     std::vector<fs::path>& Filenames);

    /**
     *  @brief    write sub catalogs
     *  @param    Filenames <vector<string> >, name of the subcatalog
     *  @param    Data <vector<double> >, data to write
     *  @param    colname <string>, column name of the respective data
     */
    void saveSubCatalogs(const std::string& filename, std::vector<double>& data,
                         const std::string &colname);

    /**
     *  @brief    returns indices of selected galaxies according to z bin
     *  @param    indices <vector<int> >, indices of selected galaxies
     *  @param    t_zmin <double>, min z value
     *  @param    t_zmax <double>, max z value
     *  @param    izbin <size_t> bin number
     */
    void getIndices(std::vector<int>& indices, double t_zmin, double t_zmax,
                    int izbin);

private:

    /**
     *  @brief parameters
     */
    const GenericParam& m_param;

    /**
     *  @brief CatalogData
     */
    CatalogData& m_cat;

    /**
     *  @brief <zMin>, Min Redshift value based on input parameter file
     *         <zMax>, Max Redshift value based on input parameter file
     */
    double m_zMin, m_zMax;

    /**
     *  @brief <m_catType>, origin of input catalog (helps to re-write sub-catalogues)
     */
    std::string m_catType;

};
// End of SplitCatalog class

/**
 *  @brief    Template to split a vector
 *  @param    vec <vector<T> >, vector to split
 *  @param    n <uint64_t>, number of sub-vectors
 *  @param    i <size_t> index of subvector to get
 *  @return   splitted vector
 */
template<typename T>
std::vector<T> split(const std::vector<T>& vec, uint64_t n, size_t i)
{
    std::vector<T> out_vec;
    uint64_t quotient = vec.size() / n;
    uint64_t reminder = vec.size() % n;

    auto start = i * quotient;
    auto stop = i * quotient + quotient;

    if (i < reminder)
    {
        stop += 1 + i;
        if (i > 0)
        {
            start += i;
        }
    }
    else
    {
        start += reminder;
        stop += reminder;
    }

    auto start_itr = std::next(vec.cbegin(), start);
    auto end_itr = std::next(vec.cbegin(), stop);
    std::copy(start_itr, end_itr, back_inserter(out_vec));

    return out_vec;
}

/**
 *  @brief  Extract elements from vector at indices specified by range
 */
template<typename ForwardIt, typename IndicesForwardIt>
inline std::vector<typename std::iterator_traits<ForwardIt>::value_type> extract_at(
        ForwardIt first, IndicesForwardIt indices_first,
        IndicesForwardIt indices_last)
{
    typedef std::vector<typename std::iterator_traits<ForwardIt>::value_type> vector_type;
    vector_type extracted;
    extracted.reserve(
            static_cast<typename vector_type::size_type>(std::distance(
                    indices_first, indices_last)));
    for (; indices_first != indices_last; ++indices_first)
        extracted.push_back(*(first + *indices_first));
    return extracted;
}

/**
 *  @brief  Extract elements from collection specified by collection of indices
 */
template<typename TVector, typename TIndicesVector>
inline TVector extract_at(const TVector& data, const TIndicesVector& indices)
{
    return extract_at(data.begin(), indices.begin(), indices.end());
}

}  // namespace LE3_2D_MASS_WL_UTILITIES

#endif
