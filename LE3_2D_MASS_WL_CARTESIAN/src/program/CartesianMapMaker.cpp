/**
 * @file src/program/LE3_2D_MASS_WL_CartesianMapMaker.cpp
 * @date 10/22/19
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

#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergencePatch.h"

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/program_options.hpp>
#include "ElementsKernel/ProgramHeaders.h"

#include "LE3_2D_MASS_WL_UTILITIES/PatchDef.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianAlgoKS.h"
#include "LE3_2D_MASS_WL_UTILITIES/NoisyCatalogData.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

using boost::program_options::options_description;
using boost::program_options::variable_value;
using std::chrono::duration;
using std::chrono::high_resolution_clock;
using std::string;
using std::vector;

using namespace LE3_2D_MASS_WL_CARTESIAN;
using namespace LE3_2D_MASS_WL_UTILITIES;

using namespace dpd::le3::wl::twodmass::inp::paramsconvergencepatch;

namespace po = boost::program_options;
namespace fs = boost::filesystem;

static Elements::Logging logger = Elements::Logging::getLogger(
        "CartesianMapMaker");

class CartesianMapMaker: public Elements::Program
{

public:

    OptionsDescription defineSpecificProgramOptions() override
    {
        OptionsDescription config_options{"Program options"};
        auto  add = config_options.add_options();

        bool flag = false;

        // Add the specific program options
        add("int-option", po::value<int>()->default_value(int{111}), "An example int option");

        return config_options;
    }

    Elements::ExitCode mainMethod(std::map<std::string, variable_value>& args)
            override
    {

        //Start time and intro
        auto start_time = std::chrono::system_clock::now();
        logger.info("#");
        logger.info("# Entering mainMethod()");
        logger.info("#");

        ////////////////////////////////////////////////////////////////////////
        // get out time running ID for filename part
        ////////////////////////////////////////////////////////////////////////
        logger.info() << "Using " << omp_get_max_threads() << " thread(s) max.";

        ////////////////////////////////////////////////////////////////////////
        // Get the workdir, the data_dir and manage Input output
        ////////////////////////////////////////////////////////////////////////
        fs::path workdir { args["workdir"].as<string>() };
        fs::path inCatalogShear { args["input_ShearCatalog"].as<string>() };
        fs::path inCatalogCluster { args["input_ClusterCatalog"].as<string>() };
        fs::path paramFile { args["paramFile"].as<string>() };
        fs::path outShearMapFilename { args["outShearMap"].as<std::string>() };
        fs::path outShearMapJson { args["outShearMapList"].as<string>() };
        fs::path datadir { workdir / "data" };

        ////////////////////////////////////////////////////////////////////////
        // get Shear catalog Data
        ////////////////////////////////////////////////////////////////////////
        CatalogData cat;
        cat.getCatalogData(workdir, inCatalogShear);

        ////////////////////////////////////////////////////////////////////////
        // check parameter file exists
        ////////////////////////////////////////////////////////////////////////
        if (false == fs::is_regular(workdir / paramFile))
        {
            logger.info()
                    << "Parameter file not found or not in XML format ....";
            return Elements::ExitCode::NOINPUT;
        }

        ////////////////////////////////////////////////////////////////////////
        // read parameter file (either convergencepatch or cluster)
        ////////////////////////////////////////////////////////////////////////
        CartesianParam params;

        // in case of patchToSphere
        double catRamin, catRamax, catDecmin, catDecmax;
        vecMinMax(cat["ra"], catRamin, catRamax);
        vecMinMax(cat["dec"], catDecmin, catDecmax);

        // in case of cluster
        CatalogData clusterCat;
        if (inCatalogCluster.string().empty() == false)
        {
            clusterCat.getCatalogData(workdir, inCatalogCluster);
            logger.info() << "cluster cat size: " << clusterCat.getNentries();
        }

        readParameterFile((workdir / paramFile), params, cat,
                           catRamin, catRamax, catDecmin, catDecmax,
                           clusterCat);

        if (outShearMapJson.string().empty())
        {
            outShearMapJson = fs::path(
                    "ShearMapsList_" + getDateTimeString() + ".json");
        }

        ////////////////////////////////////////////////////////////////////////
        // Creating patch of shear map from input catalog
        ////////////////////////////////////////////////////////////////////////
        std::ofstream outfile;
        outfile.open((workdir / outShearMapJson).string(), std::ios_base::app);
        outfile << "[";

        // iterate over each patch defined in params
        CartesianAlgoKS CartesainAlgo(params);
        for (size_t it = 0; it < params.getNPatches(); it++)
        {
            auto& patch = params.getPatches()[it];

            logger.info() << "ra_min: " << patch.getRaMin()
                          << "ra_max: " << patch.getRaMax();
            logger.info() << "dec_min: " << patch.getDecMin()
                          << "dec_max: " << patch.getDecMax();
            //logger.info() << "z_min: " << patch.getZMin()
            //              << "z_max: " << patch.getZMax();

            logger.info("creating shear map");
            if ((outShearMapFilename.string().empty()))
            {
                outShearMapFilename = fs::path(
                        "ShearMap_0" + std::to_string(it) + "_"
                                + getDateTimeString() + ".fits");
            }

            CartesainAlgo.extractShearMap(datadir / outShearMapFilename, cat, patch);
            outfile << outShearMapFilename.filename();
            if (it < params.getNPatches() - 1)
            {
                outfile << ",";
            }
            outShearMapFilename.clear();

            // Get associated SNR maps
            logger.info("creating associated snr shear Maps");
            if (params.getNResamples() > 0)
            {
                outfile << ",";
                for (int iter = 0; iter < params.getNResamples(); ++iter)
                {
                    NoisyCatalogData randomise;
                    CatalogData randCat = randomise.getNoisyCatalog(cat);
                    outShearMapFilename = fs::path(
                            "ShearMap_0" + std::to_string(it) + "_NReSample_0"
                                    + std::to_string(iter) + "_"
                                    + getDateTimeString() + ".fits");
                    CartesainAlgo.extractShearMap(datadir / outShearMapFilename, randCat,
                            patch);

                    ////////////////////////////////////////////////////////////
                    // Generate the json/txt product
                    ////////////////////////////////////////////////////////////
                    outfile << outShearMapFilename.filename();
                    if (iter < params.getNResamples() - 1)
                    {
                        outfile << ",";
                    }
                    outShearMapFilename.clear();

                } //end of if
            } // end of for loop for Nsample
        } //end of for loop

        outfile << "]";
        outfile.close();

        ////////////////////////////////////////////////////////////////////////
        //End of patch extraction
        ////////////////////////////////////////////////////////////////////////
        auto end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_time - start_time;
        logger.info() << "elapsed time: " << elapsed_seconds.count() << "s";
        logger.info("#");
        logger.info("# Exiting mainMethod()");
        logger.info("#");

        return Elements::ExitCode::OK;
    }
};

MAIN_FOR(CartesianMapMaker)
