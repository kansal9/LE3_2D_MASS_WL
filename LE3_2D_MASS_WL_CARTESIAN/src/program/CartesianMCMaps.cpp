/**
 * @file src/program/LE3_2D_MASS_WL_CartesianMCMaps.cpp
 * @date 10/25/19
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
#include "omp.h"

// Datamodel for INPUT products
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergencePatch.h"

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/program_options.hpp>
#include "ElementsKernel/ProgramHeaders.h"

#include "LE3_2D_MASS_WL_UTILITIES/PatchDef.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GetCartesianMCMaps.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianAlgoKS.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

using boost::program_options::options_description;
using boost::program_options::variable_value;
using std::chrono::duration;
using std::chrono::high_resolution_clock;
using std::string;
using std::vector;

using namespace dpd::le3::wl::twodmass::inp::paramsconvergencepatch;

using namespace LE3_2D_MASS_WL_CARTESIAN;
using namespace LE3_2D_MASS_WL_UTILITIES;

namespace fs = boost::filesystem;
namespace po = boost::program_options;

static Elements::Logging logger = Elements::Logging::getLogger(
        "CartesianMCMaps");

class CartesianMCMaps: public Elements::Program
{

public:

    options_description defineSpecificProgramOptions() override
    {

        options_description options
        { };

        // input working directory
        options.add_options()("workdir", po::value<string>()->default_value(""),
                "working directory path");

        // input log directory
        options.add_options()("logdir", po::value<string>()->default_value(""),
                "log directory path");

        // input shear catalog file
        options.add_options()("inputShearCatalog",
                po::value<string>()->default_value(""), "input shear catalog");

        // input parameter file
        options.add_options()("PatchparameterFile",
                po::value<string>()->default_value(""), "input parameter file");

        // output product (list of MC Maps)
        options.add_options()("MCConvergenceMaps",
                po::value<string>()->default_value(""),
                "MC convergence Maps Name in txt/jason file");

        // output product (list of MC Maps)
        options.add_options()("MCConvergenceMapsXML",
                po::value<string>()->default_value("outMCConvMapXML.xml"),
                "MC convergence Maps Name in XML file");

        return options;
    }

    Elements::ExitCode mainMethod(std::map<std::string, variable_value>& args)
            override
    {
        //Start time and intro
        auto start_MCMainTime = std::chrono::system_clock::now();
        logger.info("# Entering mainMethod()");
        logger.info() << "Using " << omp_get_max_threads() << " thread(s) max.";

        ////////////////////////////////////////////////////////////////////////
        // Get the workdir, the data_dir and manage Input output
        ////////////////////////////////////////////////////////////////////////
        fs::path workdir { args["workdir"].as<string>() };
        fs::path datadir { workdir / "data" };
        fs::path inFilename { args["inputShearCatalog"].as<string>() };
        fs::path inPatchParamFile
        { args["PatchparameterFile"].as<std::string>() };
        fs::path logPath(args["logdir"].as<string>());
        std::vector<fs::path> MCConvergenceMaps;
        std::vector<fs::path> MCShearMaps;
        fs::path outConv_filenames = args["MCConvergenceMaps"].as<string>();
        auto out_xml_product_file_name = args["MCConvergenceMapsXML"].as<
                std::string>();
        if ((outConv_filenames.string()).empty())
        {
            outConv_filenames = fs::path(
                    "ConvergenceMapsList_" + getDateTimeString() + ".json");
        }

        ////////////////////////////////////////////////////////////////////////
        // check parameter file exists
        ////////////////////////////////////////////////////////////////////////
        if (false == fs::is_regular(workdir / inPatchParamFile))
        {
            logger.info() << "Parameter file not found or not in XML format...";
            return Elements::ExitCode::NOINPUT;
        }

        ////////////////////////////////////////////////////////////////////////
        // Reading input catalog
        ////////////////////////////////////////////////////////////////////////
        CatalogData cat;
        logger.info("Reading Input Catalog");

        ////////////////////////////////////////////////////////////////////////
        // check whether parameter file is of XML format or not
        ////////////////////////////////////////////////////////////////////////
        CartesianParam params;
        readParameterFile((workdir / inPatchParamFile), params, cat);

        ////////////////////////////////////////////////////////////////////////
        // Variable to save input catalog name
        ////////////////////////////////////////////////////////////////////////
        std::string inCatalogFile;

        ////////////////////////////////////////////////////////////////////////
        // Loading shear catalog
        ////////////////////////////////////////////////////////////////////////
        if (inFilename.string().empty() == false)
        {
            cat.getCatalogData(workdir, inFilename);
        }
        else
        {
            logger.info() << "ERROR: There is no Input Shear Catalogue...";
            return Elements::ExitCode::NOINPUT;
        }

        ////////////////////////////////////////////////////////////////////////
        // get Denoised Shear Map
        ////////////////////////////////////////////////////////////////////////
        logger.info() << "get Denoised Shear Map";

        double ramin, decmin;
        GetCartesianMCMaps MCMaps(cat, params);

        for (size_t it = 0; it < params.getNPatches(); it++)
        {
            auto& patch = params.getPatches()[it];
            ramin = patch.getRaMin();
            decmin = patch.getDecMin();

            ShearMap deNoisedShearMap;
            MCMaps.getDeNoisedShearMap(deNoisedShearMap);
            auto deNoisedShearMapPath = datadir
                    / fs::path(
                           "DenoisedShearMap_" + getDateTimeString() + ".fits");
            logger.info() << "Will write denoisedShearMap to"
                    << deNoisedShearMapPath;
            deNoisedShearMap.writeMap(deNoisedShearMapPath, params);

            ////////////////////////////////////////////////////////////////////
            // get noised N Shear Maps
            ////////////////////////////////////////////////////////////////////
            logger.info() << "get N noised shear maps";

            std::vector<ShearMap> NoisedShearMapList;
            MCMaps.getNoisedShearMaps(NoisedShearMapList);

            for (size_t i = 0; i < NoisedShearMapList.size(); i++)
            {
                fs::path shearMapName = (datadir
                        / fs::path(
                                "EUC_LE3_WL_NoisedShearMap_" + std::to_string(i)
                                        + "_" + getDateTimeString() + ".fits"));
                MCShearMaps.push_back(shearMapName);
            }

            for (unsigned int i = 0; i < NoisedShearMapList.size(); i++)
            {
                logger.info() << i;
                MCMaps.performAddition(deNoisedShearMap, NoisedShearMapList[i],
                        MCShearMaps[i].native());
            }

            ////////////////////////////////////////////////////////////////////
            // perform mass mapping on N shear Maps
            ////////////////////////////////////////////////////////////////////
            logger.info("entering Mass mapping()");
            for (size_t i = 0; i < MCShearMaps.size(); i++)
            {
                fs::path convMapFilename = fs::path(
                        "EUC_LE3_WL_MCConvergenceMap_" + std::to_string(i) + "_"
                                + getDateTimeString() + ".fits");
                fs::path convMapFilepath = datadir / convMapFilename;
                ShearMap shearMap(MCShearMaps[i]);
                ConvergenceMap convergenceMap(shearMap);
                shearMap.getConvMap(convergenceMap);
                convergenceMap.writeMap(convMapFilepath, params);
                MCConvergenceMaps.push_back(convMapFilepath);
            }

            ////////////////////////////////////////////////////////////////////
            // Generate the json/txt product
            ////////////////////////////////////////////////////////////////////
            fillJsonOutput(workdir, outConv_filenames, MCConvergenceMaps);
        }

        ////////////////////////////////////////////////////////////////////////
        //End of MC Main
        ////////////////////////////////////////////////////////////////////////
        auto end_MCMainTime = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_MCMainTime
                - start_MCMainTime;
        logger.info() << "elapsed time: " << elapsed_seconds.count() << "s";
        logger.info("#");
        logger.info("# Exiting mainMethod()");
        logger.info("#");
        return Elements::ExitCode::OK;
    }
};

MAIN_FOR(CartesianMCMaps)
