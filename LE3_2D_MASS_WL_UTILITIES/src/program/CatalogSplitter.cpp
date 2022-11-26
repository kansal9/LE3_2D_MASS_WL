/**
 * @file src/program/LE3_2D_MASS_WL_CatalogSplitter.cpp
 * @date 10/23/19
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

#include <map>
#include <string>
#include <chrono>
#include <ctime>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include "omp.h"

// Datamodel for INPUT products
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergencePatch.h"

#include <boost/program_options.hpp>
#include "ElementsKernel/ProgramHeaders.h"

#include "LE3_2D_MASS_WL_UTILITIES/SplitCatalogRedshift.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

using boost::program_options::options_description;
using boost::program_options::variable_value;
using std::chrono::duration;
using std::chrono::high_resolution_clock;
using std::string;
using std::vector;

using namespace LE3_2D_MASS_WL_UTILITIES;
using LE3_2D_MASS_WL_CARTESIAN::CartesianParam;
using LE3_2D_MASS_WL_SPHERICAL::SphericalParam;

static Elements::Logging logger =
                                Elements::Logging::getLogger("CatalogSplitter");

/**
 *  @brief    Template to sort indices of a vector
 *  @param    v <vector<T> >, vector to sort
 *  @return   vector of sorted indices
 */
template<typename T>
std::vector<size_t> getSortedIndex(const std::vector<T> &v)
{
    // initialize original index
    std::vector<size_t> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i)
    {
        idx[i] = i;
    }
    // sort indexes based by comparing values in v using std::stable_sort
    // instead of std::sort to avoid unnecessary index re-orderings when v
    // contains elements of equal values
    sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2)
    {
        return v[i1] < v[i2];
    });
    return idx;
}

void sortCatalogData(CatalogData& cat)
{
    // init temporary catalog
    CatalogData temp(cat);

    // get sorted index according to reference column for sort
    auto v = cat["z"].container();
    auto index = getSortedIndex(v);

    // apply sort to catalog
    for(auto const& col: cat.getColumnProxyNames())
    {
        for(long j=0; j<cat.getNentries(); j++)
        {
            cat[col](j) = temp[col](index[j]);
        }
    }
}


class CatalogSplitter: public Elements::Program
{

public:

    options_description defineSpecificProgramOptions() override
    {

        options_description options
        { };
        // input working directory
        options.add_options()("workdir", po::value<string>()->default_value(""),
                "Working Directory Path");

        // input shear catalog file
        options.add_options()("inputCatalog",
                po::value<string>()->default_value(""), "input shear Catalog");

        // input parameter file
        options.add_options()("paramFile",
                po::value<string>()->default_value(""),
                "Input Parameter File to split input catalog");

        // output product (list of sub-catalogs)
        options.add_options()
        ("Sub_Catalogs", po::value<string>()->default_value("SubCatalogs.json"),
                "Sub-Catalog Name in txt/jason file");

        return options;
    }

    Elements::ExitCode mainMethod(std::map<std::string, variable_value>& args)
    override
    {
        //Start time and intro
        auto start_SplitterTime = std::chrono::system_clock::now();
        logger.info("#");
        logger.info("# Entering mainMethod()");
        logger.info("#");

        ////////////////////////////////////////////////////////////////////////
        // get out time running ID for filename part
        ////////////////////////////////////////////////////////////////////////
        // logger.info() << "Using " << omp_get_max_threads() << " thread(s) max.";

        ////////////////////////////////////////////////////////////////////////
        // Get the workdir, the datadir and manage Input output
        ////////////////////////////////////////////////////////////////////////
        fs::path workdir { args["workdir"].as<string>() };
        fs::path datadir { workdir / "data" };
        fs::path inputCatalog { args["inputCatalog"].as<string>() };
        fs::path parameterFile { args["paramFile"].as<std::string>() };
        std::vector<fs::path> SubCatFilenames;

        ////////////////////////////////////////////////////////////////////////
        // check parameter file exists
        ////////////////////////////////////////////////////////////////////////
        if (false == fs::is_regular(workdir / parameterFile))
        {
            logger.info()
                    << "Parameter file not found or not in XML format ....";
            return Elements::ExitCode::NOINPUT;
        }

        ////////////////////////////////////////////////////////////////////////
        // Reading input catalog
        ////////////////////////////////////////////////////////////////////////
        CatalogData cat;
        cat.getCatalogData(workdir, inputCatalog);

        ////////////////////////////////////////////////////////////////////////
        // check whether parameter file is of XML format
        ////////////////////////////////////////////////////////////////////////
        GenericParam params;

        std::string xmlProductType = getXmlProductType(parameterFile);
        if(xmlProductType == "DpdTwoDMassParamsConvergencePatch"
        or xmlProductType == "DpdTwoDMassParamsConvergencePatchesToSphere")
        {
            CartesianParam& cartesianParam = static_cast<CartesianParam&>(params);
            readParameterFile(parameterFile, cartesianParam, cat);
        }
        else if(xmlProductType == "DpdTwoDMassParamsConvergenceSphere")
        {
            SphericalParam& sphericalParam = static_cast<SphericalParam&>(params);
            sphericalParam.readConvSphereXMLFile(parameterFile, cat);
        }

        std::string catType = cat.getCatalogType();
        fs::path inputCatalogName = cat.getCatalogFitsName();
        if (params.getNbins() <= 1)
        {
            SubCatFilenames.push_back(inputCatalogName);
        }

        if (params.getNbins() > 1)
        {
            ////////////////////////////////////////////////////////////////////
            // If redshift is 0 for all objects in catalog and balanced bin is
            // false, then set balanced bin to true
            ////////////////////////////////////////////////////////////////////
            double zMin, zMax;
            vecMinMax(cat["z"], zMin, zMax);

            if (zMin == 0 && zMax == 0)
            {
                logger.info()
                      << "ERROR: The redshift values for all objects inside catalog is 0.";
                return Elements::ExitCode::OK;
            }

            for (int i = 0; i < params.getNbins(); i++)
            {
                fs::path subCatalogName = fs::path(
                        "EUC_LE3_WL_CatZ_subcatalogue_Z0" + std::to_string(i)
                                + "_" + getDateTimeString() + ".fits");
                SubCatFilenames.push_back(subCatalogName);
            }

            if (params.isBalancedBins())
            {
                // sort according to redshift
                sortCatalogData(cat);
            }

            ////////////////////////////////////////////////////////////////////////
            // Splitting Input catalog
            ////////////////////////////////////////////////////////////////////////
            SplitCatalogRedshift split(cat, params, catType);
            split.writeSubCatalogs(datadir, SubCatFilenames);
        }
        ////////////////////////////////////////////////////////////////////////
        // Generate the jason/txt product
        ////////////////////////////////////////////////////////////////////////
        fs::path output_filename = args["Sub_Catalogs"].as<string>();
        fillJsonOutput(workdir, output_filename, SubCatFilenames);

        ////////////////////////////////////////////////////////////////////////
        //End of Splitter
        ////////////////////////////////////////////////////////////////////////
        auto end_SplitterTime = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_SplitterTime
                - start_SplitterTime;
        logger.info() << "elapsed time: " << elapsed_seconds.count() << "s";
        logger.info("#");
        logger.info("# Exiting mainMethod()");
        logger.info("#");
        return Elements::ExitCode::OK;
    }

};

MAIN_FOR(CatalogSplitter)
