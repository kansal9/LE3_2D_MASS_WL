/**
 * @file src/lib/CatalogData.cpp
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

#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"

using namespace Euclid;
using namespace Fits;

using Elements::ExitCode;

static Elements::Logging logger = Elements::Logging::getLogger("CatalogData");

namespace LE3_2D_MASS_WL_UTILITIES {

CatalogData::CatalogData() : m_inputCatalog(), m_catType(), m_nEntries(0)
{
}

CatalogData::CatalogData(const CatalogData& cat) :
        m_inputCatalog(cat.m_inputCatalog),
        m_catType(cat.m_catType),
        m_nEntries(cat.m_nEntries), m_data(cat.m_data),
        m_colNamesHdu(cat.m_colNamesHdu), m_colNames(cat.m_colNames),
        m_sortedIndex(cat.m_sortedIndex)
{
}

void CatalogData::getCatalogData(fs::path& workdir, fs::path& inputCatalog)
{
    fs::path datadir = workdir / "data";
    // Case 1: When input catalog file is in fits format:
    // directly read this catalog (expected workdir/data)
    if (checkFileType(datadir / inputCatalog, signFITS))
    {
        m_inputCatalog = datadir / inputCatalog;
    }
    // Case 2: When input catalog file is in xml format:
    // fetch catalog filename from xml file and read it (expected in workdir)
    else if (checkFileType(workdir / inputCatalog, signXML))
    {
        m_inp.readXmlFile(workdir / inputCatalog);
        m_inputCatalog = datadir / m_inp.getCatalogFilename();
    }
    readCatalog(m_inputCatalog);
    logger.info() << "Done reading catalog: " << m_inputCatalog
                  << " (nEntries=" << getNentries() << ")";

}

void CatalogData::readCatalog(const fs::path& filename)
{
    if (!checkFileType(filename, signFITS))
    {
        logger.fatal() << "Bad catalog file: ", filename.native();
        throw ExitCode::NOINPUT;
    }
    MefFile f(filename.native(), FileMode::Read);
    // access hdu indexed 0 (not using its name)
    const BintableHdu& ext = f.access<BintableHdu>(1);
    const BintableColumns& cols = ext.columns();
    m_colNamesHdu = cols.readAllNames();
    m_nEntries = cols.readRowCount();
    setColnamesToRead();
    std::vector<ColumnKey> colKeys;
    for(size_t i = 0; i < m_colNames.size(); i++)
    {
        colKeys.push_back(ColumnKey(m_colNames[i].first));
    }
    m_data = cols.readSeq<double>(colKeys);
    m_sortedIndex = std::vector<int>(m_nEntries);
    std::iota(m_sortedIndex.begin(), m_sortedIndex.end(), 0);
}

std::set<std::string> CatalogData::getShearTypesInCat()
{
    std::set<std::string> shearTypes;
    const std::regex col_template("SHE_(.*)_G1");
    std::smatch match;
    for (const auto &name : m_colNamesHdu)
    {
        if (std::regex_match(name, match, col_template))
        {
            // First sub_match is the whole string; next is the first group
            if (match.size() == 2)
            {
                shearTypes.insert(match[1].str());
            }
        }
    }
    return shearTypes;
}

void CatalogData::setColnamesToRead()
{
    m_colNames.clear();
    std::set<std::string> shearAvailTypes = getShearTypesInCat();

    m_catType = "LENSMC";
    if(shearAvailTypes.size() == 1)
    {
        m_catType = *(shearAvailTypes.begin());
        m_colNames = getShearColNamesAndProxy(m_catType);
    }
    else if(shearAvailTypes.size() > 1)
    {
        // here we can choose a given shear type within the possible ones
        // TODO: add a parameter to choose one different from LENSMC if possible
        std::set<std::string>::iterator it = shearAvailTypes.find(m_catType);
        if(it != shearAvailTypes.end())
        {
            m_catType = *it;
        }
        // we are on a shear catalog with a type of shear given by shearType
        m_colNames = getShearColNamesAndProxy(m_catType);
    }
    else
    {
        // else its a cluster one
        m_catType = "CLUSTER";
        m_colNames = getClusterColNamesAndProxy();
    }
}

double CatalogData::operator()(const std::string& colname, int i) const
{
    int index = getColumnIndex(colname);
    if(i < 0 or i > m_nEntries)
    {
        logger.fatal() << "Index not found: " << i;
    }
    return m_data[index](m_sortedIndex[i]);
}

const VecColumn<double>& CatalogData::operator[](const std::string& colname) const
{
    int index = getColumnIndex(colname);
    return m_data[index];
}

VecColumn<double>& CatalogData::operator[](const std::string& colname)
{
    int index = getColumnIndex(colname);
    return m_data[index];
}

void CatalogData::fillTest(int N, const std::string& shearType,
                           double raMin, double raMax,
                           double decMin, double decMax,
                           double zMin, double zMax)
{
    logger.info() << "Fill for test purpose with patch bounds: "
                  << raMin << " " << raMax << " "
                  << decMin << " " << decMax << " "
                  << zMin << " " << zMax;

    m_colNames = getShearColNamesAndProxy(shearType);
    m_nEntries = N;

    m_data.resize(m_colNames.size());
    for(size_t i = 0; i < m_colNames.size(); i++)
    {
        m_data[i] = VecColumn<double, 1>(m_data[i].info(), N);
        m_data[i].rename(m_colNames[i].first);
    }
    fillVecColumnLinspace((*this)["ra"], N, raMin, raMax);
    fillVecColumnLinspace((*this)["dec"], N, decMin, decMax);
    fillVecColumnLinspace((*this)["w"], N, 1, 1);
    fillVecColumnLinspace((*this)["g1"], N, 0, 0.4);
    fillVecColumnLinspace((*this)["g2"], N, 0, 0.1);
    fillVecColumnLinspace((*this)["z"], N, zMin, zMax);
    fillVecColumnLinspace((*this)["z_corr"], N, 0, 0);
    m_sortedIndex = std::vector<int>(m_nEntries);
    std::iota(m_sortedIndex.begin(), m_sortedIndex.end(), 0);
}

long CatalogData::getNentries() const
{
    return m_nEntries;
}

size_t CatalogData::getNcols() const
{
    return m_colNames.size();
}

fs::path CatalogData::getCatalogFitsName() const
{
    return m_inputCatalog.filename();
}

std::string CatalogData::getCatalogType() const
{
    return m_catType;
}

std::vector<std::pair<std::string, std::string>>
CatalogData::getColumnNamePairs() const
{
    return m_colNames;
}

std::vector<std::string> CatalogData::getColumnProxyNames() const
{
    std::vector<std::string> columnProxyNames;
    for(auto const& col: m_colNames)
    {
        columnProxyNames.push_back(col.second);
    }
    return columnProxyNames;
}

void CatalogData::sortIndex(const std::string& key)
{
    int index = getColumnIndex(key);
    sort(m_sortedIndex.begin(), m_sortedIndex.end(), [this, index](int i, int j)
    {
        return m_data[index](i) < m_data[index](j);
    });
}

void CatalogData::getMinMax(const std::string& key, double& min, double& max) const
{
    int index = getColumnIndex(key);

    auto res = std::minmax_element(m_data[index].begin(), m_data[index].end());
    min = *res.first;
    max = *res.second;
}


int CatalogData::getColumnIndex(std::string colname) const
{
    // search for colname in proxy name (second entry of pair)
    auto it = std::find_if(m_colNames.begin(), m_colNames.end(),
        [colname](const std::pair<std::string, std::string>& entry)
        { return entry.first == colname or entry.second == colname;} );
    if(it == m_colNames.end())
    {
        logger.error() << "Column name not found (return index 0): " << colname;
        return 0;
    }
    else
    {
        return std::distance(m_colNames.begin(), it);
    }
}

long CatalogData::getNentriesBounds(const std::string& col, double lo, double hi) const
{
    int index = getColumnIndex(col);
    return std::count_if(m_data[index].begin(), m_data[index].end(),
            [lo, hi](double v) { return (lo <= v) and (v <= hi); });
}

}  // namespace LE3_2D_MASS_WL_UTILITIES

