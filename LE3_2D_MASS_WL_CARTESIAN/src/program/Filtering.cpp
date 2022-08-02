/**
 * @file src/program/LE3_2D_MASS_WL_Filtering.cpp
 * @date 04/06/21
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

#include <boost/program_options.hpp>
#include "ElementsKernel/ProgramHeaders.h"

#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergencePatch.h"

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/MatrixProcess.h"
#include "LE3_2D_MASS_WL_CARTESIAN/FilterMR.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GenericMap.h"

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>

#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <iomanip>

namespace po = boost::program_options;
using std::string;
using std::vector;

using boost::program_options::options_description;
using boost::program_options::variable_value;

using namespace dpd::le3::wl::twodmass::inp::paramsconvergencepatch;
using namespace LE3_2D_MASS_WL_UTILITIES;
using namespace LE3_2D_MASS_WL_CARTESIAN;

static Elements::Logging logger = Elements::Logging::getLogger(
        "Filtering");

class Filtering: public Elements::Program
{

public:

    options_description defineSpecificProgramOptions() override
    {

        options_description options
        { };
        // input working directory
        options.add_options()("workdir", po::value<string>()->default_value(""),
                "Work Directory");

        // input log directory: default is none
        options.add_options()("logdir", po::value<string>()->default_value(""),
                "logs Directory");

        // input catalog file
        options.add_options()("inmap", po::value<string>()->default_value(""),
                "input map");

        // input parameter file
        options.add_options()("paramFile",
                po::value<string>()->default_value(""),
                "Input Parameter File in XML");

        options.add_options()("positiveCons",
                po::value<int>()->default_value(0),
                "Positive constraint (0-> False and 1-> True)");

        options.add_options()("KillLastScale",
                po::value<int>()->default_value(0),
                "Kill last scale during the filtering (0-> False and 1-> True)");

        options.add_options()("KillIsol", po::value<int>()->default_value(0),
                "suppress isolated pixels (0-> False and 1-> True)");

        options.add_options()("FirstScale", po::value<int>()->default_value(0),
                "FirstScale used for detection");

        options.add_options()("nbIter", po::value<int>()->default_value(10),
                "number of loops in reconstruction iterative process");

        // output file
        options.add_options()("outputMap",
                po::value<string>()->default_value(""),
                "Output Map in fits format");
        return options;
    }

    Elements::ExitCode mainMethod(std::map<std::string, variable_value>& args)
    override
    {
        logger.info("#");
        logger.info("# Entering mainMethod()");
        logger.info("#");

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // get out time running ID for filename part
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        logger.info() << "Using " << omp_get_max_threads() << " thread(s) max.";

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Get the workdir, the data_dir and manage Input output
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        fs::path workdir
        { args["workdir"].as<string>() };
        fs::path Inmap
        { args["inmap"].as<string>() };
        fs::path ParameterFile
        { args["paramFile"].as<std::string>() };
        fs::path logPath(args["logdir"].as<string>());
        fs::path DenoisedMap
        { args["outputMap"].as<std::string>() };
        fs::path datadir
        { workdir / "data" };
        bool positiveCons = true;
        bool KillLastScale = true;
        bool KillIsol = true;
        int nbIter
        { args["nbIter"].as<int>() };
        int FirstScale
        { args["FirstScale"].as<int>() };

        if (FirstScale <= 0)
        {
            FirstScale = 1;
        }

        if (args["positiveCons"].as<int>() == 0)
        {
            positiveCons = false;
        }

        if (args["KillLastScale"].as<int>() == 0)
        {
            KillLastScale = false;
        }

        if (args["KillIsol"].as<int>() == 0)
        {
            KillIsol = false;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // check parameter file exists
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        if (false == fs::is_regular(workdir / ParameterFile))
        {
            logger.info()
                    << "Parameter file not found or not in XML format ....";
            return Elements::ExitCode::NOINPUT;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // read parameter File
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        CartesianParam params;
        CatalogData dummy;
        readParameterFile((workdir / ParameterFile), params, dummy);
        logger.info() << "FDR_val = " << params.getThresholdFdr();

        LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap convergenceMap(
                (workdir / Inmap).native());

        FilterMR filter(convergenceMap, positiveCons, KillLastScale, KillIsol,
                nbIter, FirstScale, params);

        filter.performFiltering(convergenceMap);

        std::string timeNamePart = getDateTimeString();
        if (true == DenoisedMap.string().empty())
            DenoisedMap = fs::path(
                    "EUC_LE3_WL_DenoisedMap_" + timeNamePart + ".fits");
        convergenceMap.writeMap((workdir / DenoisedMap).native(), params);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        //End of Filtering
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        logger.info("#");
        logger.info("# Exiting mainMethod()");
        logger.info("#");

        return Elements::ExitCode::OK;
    }

};

MAIN_FOR(Filtering)
