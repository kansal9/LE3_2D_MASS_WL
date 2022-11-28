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
#include <regex>
#include <boost/filesystem.hpp>

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"

#include "ElementsKernel/Logging.h"

#include "EleFits/MefFile.h"

namespace fs = boost::filesystem;

using LE3_2D_MASS_WL_UTILITIES::DmInput;

namespace LE3_2D_MASS_WL_UTILITIES {

/**
 * @class CatalogData
 * @brief
 *
 */
class CatalogData
{

public:

    /**
     * @brief Destructor
     */
    virtual ~CatalogData() = default;

    /**
     @brief Constructor
     */
    CatalogData();

    /**
     @brief Copy constructor
     */
    CatalogData(const CatalogData& cat);

    /**
     * @brief         check catalog format (fits/xml) and read data
     *                of that catalog
     * @param[in]     <workdir> work directory
     * @param[in]     <inputCatalog> input catalog name
     */
    void getCatalogData(fs::path& workdir, fs::path& inputCatalog);

    /**
     * @brief         method to read the catalog
     * @param[in]     <filename> path of the fits catalog
     */
    void readCatalog(const fs::path& filename);

    /**
     * @brief    element access (read only)
     * @param[in] string colname : proxy name of the column
     * @param[in] int i : index in the column (through sorted index)
     */
    double operator()(const std::string& colname, int i) const;

    /**
     * @brief    element access (read only)
     * @param[in] string colname : proxy name of the column
     */
    const VecColumn<double>& operator[](const std::string& colname) const;

    /**
     * @brief    element access (read/write)
     * @param[in] string colname : proxy name of the column
     */
    VecColumn<double>& operator[](const std::string& colname);

    /**
     * @brief    get available shear types in catalog
     */
    std::set<std::string> getShearTypesInCat();

    /**
     * @brief    set the column names to read
     */
    void setColnamesToRead();

    /**
     * @brief    fill with test data
     */
    void fillTest(int N, const std::string& shearType="LENSMC",
                  double raMin=0, double raMax=360,
                  double decMin=0., double decMax=180,
                  double zMin=0, double zMax=9);

    /**
     * @brief    get number of entries in the catalog
     */
    long getNentries() const;

    /**
     * @brief    get number of entries in the catalog
     */
    long getNentriesBounds(const std::string& col, double lo, double hi) const;

    /**
     * @brief    get number of columns in the catalog
     */
    size_t getNcols() const;

    /**
     * @brief    get vector of column proxy names
     */
    std::vector<std::string> getColumnProxyNames() const;

    /**
     * @brief    get catalog fits name
     */
    fs::path getCatalogFitsName() const;

    /**
     * @brief    get catalog type
     */
    std::string getCatalogType() const;

    /**
      * @brief    get column names
      */
    std::vector<std::pair<std::string, std::string>> getColumnNamePairs() const;

    /**
     * @brief    sort index according of column key (ascending)
     */
    void sortIndex(const std::string& key);

    /**
     * @brief    get extremas for one given column
     */
    void getMinMax(const std::string& key, double& min, double& max) const;

    /**
     * @brief    get colupn index for column proxy name
     */
    int getColumnIndex(std::string colname) const;


private:
    /**
     * @brief    path to the fits catalog
     */
    fs::path m_inputCatalog;

    /**
     * @brief    catalog type
     */
    std::string m_catType;

    /**
     * @brief    data model input reader
     */
    DmInput m_inp;

    /**
     * @brief    number of entries
     */
    long m_nEntries;

    /**
     * @brief    data from requested columns
     */
    std::vector<VecColumn<double>> m_data;

    /**
     * @brief    all column names in the bintable hdu
     */
    std::vector<std::string> m_colNamesHdu;

    /**
     * @brief    selected column names in the bintable hdu and proxy
     */
    std::vector<std::pair<std::string, std::string>> m_colNames;

    /**
     * @brief    vector of sorted index
     */
    std::vector<int> m_sortedIndex;

}; // End of CatalogData class

}  // namespace LE3_2D_MASS_WL_UTILITIES

#endif
