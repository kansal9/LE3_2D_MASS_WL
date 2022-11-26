/**
 * @file src/lib/SplitCatalog.cpp
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

#include "LE3_2D_MASS_WL_UTILITIES/SplitCatalogRedshift.h"

static Elements::Logging logger = Elements::Logging::getLogger(
        "SplitCatalog");

namespace LE3_2D_MASS_WL_UTILITIES
{

SplitCatalogRedshift::SplitCatalogRedshift(CatalogData& cat, const GenericParam& param,
        std::string& catType) : m_param(param), m_cat(cat), m_catType(catType)
{
    double redshiftMin, redshiftMax;
    vecMinMax(m_cat["z"], redshiftMin, redshiftMax);
    double zMin = m_param.getZMin(0);
    double zMax = m_param.getZMax(0);
    if (zMin < redshiftMin)
    {
        m_zMin = redshiftMin;
    }
    else
    {
        m_zMin = zMin;
    }
    if (zMax > redshiftMax)
    {
        m_zMax = redshiftMax;
    }
    else
    {
        m_zMax = zMax;
    }
    logger.info() << " Catalog min max redshift " << m_zMin << " " << m_zMax;
}


void SplitCatalogRedshift::getEqualBinSplittedCatalogs(fs::path& datadir,
        std::vector<fs::path>& Filenames)
{
    for (int iter = 0; iter < m_param.getNbins(); iter++)
    {
        const auto filename = datadir.string()+ "/" +Filenames[iter].string();
        for (auto const& entry : m_cat.getColumnNamePairs())
        {
            std::vector<double> v_col = m_cat[entry.first].container();
            std::vector<double> temp = split(v_col, m_param.getNbins(), iter);
            saveSubCatalogs(filename, temp, entry.first);
        }
    }
}

void SplitCatalogRedshift::getSplittedCatalogs(fs::path& datadir,
        std::vector<fs::path>& Filenames)
{
    double range = (m_zMax - m_zMin) / m_param.getNbins();
    double t_zmin = m_zMin;
    double t_zmax = m_zMin + range;

    for (int it = 0; it < m_param.getNbins(); it++)
    {
        const auto filename = datadir.string() + "/" + Filenames[it].string();

        logger.info() << "Redshift Range for subcatalogue " << int(it + 1)
                        << " is z_min: " << t_zmin << "\t"
                        << "z_max: " << t_zmax;
        if (t_zmax > m_zMax)
        {
            logger.info() << " Redshift values are greater than input Zmax ";
            break;
        }

        std::vector<int> indices;
        getIndices(indices, t_zmin, t_zmax, it);
        indices.shrink_to_fit();

        logger.info() << "Length indices " << indices.size();

        for (auto const& entry : m_cat.getColumnNamePairs())
        {
            std::vector<double> v_col = m_cat[entry.first].container();
            std::vector<double> temp = extract_at(v_col, indices);
            saveSubCatalogs(filename, temp, entry.first);
        }

        logger.info() << "number of created subcatalogs: " << int(it + 1);
        t_zmin = t_zmin + range;
        t_zmax = t_zmax + range;
    }
}

void SplitCatalogRedshift::getIndices(std::vector<int>& indices, double t_zmin,
        double t_zmax, int izbin)
{
    long nGal = m_cat.getNentries();
    logger.info() << "Number of galaxies " << nGal;

    for (long i = 0; i < nGal; i++)
    {
        if (izbin < m_param.getNbins() - 1)
        {
            if (m_cat["z"](i) >= t_zmin && m_cat["z"](i) < t_zmax)
            {
                indices.push_back(i);
            }
        }
        else
        {
            if (m_cat["z"](i) >= t_zmin && m_cat["z"](i) <= t_zmax)
            {
                indices.push_back(i);
            }
        }

    }
}

void SplitCatalogRedshift::writeSubCatalogs(fs::path& datadir,
        std::vector<fs::path>& Filenames)
{
    logger.info() << "number of redshift bins are:" << m_param.getNbins();
    logger.info() << "Balanced bins:" << m_param.isBalancedBins();

    if (m_zMin < 0.0)
    {
        logger.info()
                << "Warning: Redshift Value is less than ZERO. The minimum redshift value is: "
                << m_zMin;
    }
    if (m_zMax > 10.0)
    {
        logger.info()
                << "Warning: Redshift Value is greater than TEN. The maximum redshift value is: "
                << m_zMax;
    }

    if (m_param.isBalancedBins())
    {
        logger.info() << "entering balanced bins split";
        getEqualBinSplittedCatalogs(datadir, Filenames);
    }
    else
    {
        logger.info() << "entering split";
        getSplittedCatalogs(datadir, Filenames);
    }
}

void SplitCatalogRedshift::saveSubCatalogs(const std::string& filename,
        std::vector<double>& data, const std::string &colname)
{
    VecColumn<double> col = generateColumn(colname, data);
    if (access(filename.c_str(), F_OK) != -1)
    {
        MefFile outfile(filename, FileMode::Edit);
        long nbHdu = outfile.hduCount();
        const auto &ext = outfile.access<BintableHdu>(nbHdu - 1);
        const auto &btCol = ext.columns();
        btCol.init(col.info());
        btCol.write(col);
    }
    else
    {
        MefFile outfile(filename, FileMode::Overwrite);
        outfile.assignBintableExt("", col);
    }
}

}  // namespace LE3_2D_MASS_WL_UTILITIES
