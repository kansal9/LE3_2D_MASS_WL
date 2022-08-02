/**
 * @file src/program/LE3_2D_MASS_WL_CartesianMassMapping.cpp
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
#include <vector>
#include <string>
#include <chrono>
#include <ctime>
#include "omp.h"

#include "ST_DataModelBindings/dpd/le3/wl/twodmass/out/euc-test-le3-wl-twodmass-ConvergencePatch.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergencePatch.h"

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/program_options.hpp>
#include "ElementsKernel/ProgramHeaders.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianAlgoKS.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmOutput.h"

using boost::program_options::options_description;
using boost::program_options::variable_value;
using std::chrono::duration;
using std::chrono::high_resolution_clock;
using std::string;

using namespace LE3_2D_MASS_WL_CARTESIAN;
using namespace LE3_2D_MASS_WL_UTILITIES;

using namespace dpd::le3::wl::twodmass::inp::paramsconvergencepatch;
using namespace dpd::le3::wl::twodmass::out::convergencepatch;

namespace po = boost::program_options;
namespace fs = boost::filesystem;

static Elements::Logging logger = Elements::Logging::getLogger(
        "CartesianMassMapping");

class CartesianMassMapping: public Elements::Program
{
public:

    options_description defineSpecificProgramOptions() override
    {
        options_description options{ };

        // input working directory
        options.add_options()("workdir", po::value<string>()->default_value(""),
                "working directory");

        // input log directory: default is none
        options.add_options()("logdir", po::value<string>()->default_value(""),
                "log directory");

        // input parameter file
        options.add_options()("paramFile",
                po::value<string>()->default_value(""),
                "input parameter file in xml format");

        // input shear map file
        options.add_options()("inputShearMap",
                po::value<string>()->default_value(""),
                "input shear maps in json/fits format");

        // output products (Convergence Map)
        options.add_options()("outConvMap",
                po::value<string>()->default_value(""),
                "output convergence map in fits format");
        options.add_options()("outConvMaps",
                po::value<string>()->default_value(""),
                "output list of convergence maps in json format");
        options.add_options()("outConvMapXML",
                po::value<string>()->default_value("outConvMapXML.xml"),
                "output convergence Map in xml format");

        return options;
    }

    Elements::ExitCode mainMethod(std::map<std::string, variable_value>& args)
            override
    {
        //Start time and intro
        auto start_MappingTime = std::chrono::system_clock::now();
        logger.info("# Entering mainMethod()");
        logger.info() << "Using " << omp_get_max_threads() << " thread(s) max.";

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Get the workdir, the data_dir and manage Input output
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        fs::path workdir
        { args["workdir"].as<string>() };
        fs::path datadir
        { workdir / "data" };
        fs::path paramFile
        { args["paramFile"].as<std::string>() };
        fs::path logPath(args["logdir"].as<string>());
        fs::path inShearMap
        { args["inputShearMap"].as<string>() };  // either fits or json
        fs::path outConvMap
        { args["outConvMap"].as<std::string>() };
        fs::path outConvMapsJson
        { args["outConvMaps"].as<std::string>() };
        auto out_xml_product_file_name =
                args["outConvMapXML"].as<std::string>();

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // check parameter file exists
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        if (false == fs::is_regular(workdir / paramFile))
        {
            logger.info()
                    << "Parameter file not found or not in XML format ....";
            return Elements::ExitCode::NOINPUT;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // check whether parameter file is of XML format or not
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        CartesianParam params;
        CatalogData dummy;
        readParameterFile((workdir / paramFile), params, dummy);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Object to perform Cartesian KS algorithm
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        CartesianAlgoKS CartesainAlgo(params);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // set FITS output filenames and perform KS Mass Mapping
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        // set output convergence maps list in json default name if empty
        if ((outConvMapsJson.string()).empty())
        {
            outConvMapsJson = fs::path(
                    "ConvergenceMapsList_" + getDateTimeString() + ".json");
        }

        // perform mass mapping
        CartesainAlgo.performReducedShear(inShearMap, workdir, outConvMapsJson);

        // read filenames in json file outConvMapsJson
        std::vector<fs::path> convMapsJson = readFilenamesInJson(
                outConvMapsJson);

        // if more than one entry in outConvMapsJson, compute SNR
        if (convMapsJson.size() > 1)
        {
            CartesainAlgo.getSNRMap(workdir, outConvMapsJson, outConvMap);

            // add SNR map to convMapsJson
            convMapsJson.push_back(outConvMap);
        }

        // update json
        fillJsonOutput(workdir, outConvMapsJson, convMapsJson);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        //Write XML files
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        fs::path outputXML
        { workdir / out_xml_product_file_name };
        fs::path inputjson
        { workdir / outConvMapsJson };
        CartesainAlgo.writeXMLfile(inputjson, outputXML);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        //End of Mass Mapping
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        auto end_MappingTime = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_MappingTime
                - start_MappingTime;
        logger.info() << "elapsed time: " << elapsed_seconds.count() << "s";
        logger.info("#");
        logger.info("# Exiting mainMethod()");
        logger.info("#");

        return Elements::ExitCode::OK;
    }
};

MAIN_FOR(CartesianMassMapping)
