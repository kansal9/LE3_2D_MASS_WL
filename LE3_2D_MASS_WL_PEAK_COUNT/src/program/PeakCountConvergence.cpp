/**
 * @file src/program/LE3_2D_MASS_WL_PeakCountConvergence.cpp
 * @date 12/03/19
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

#include "LE3_2D_MASS_WL_PEAK_COUNT/PeakParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ConvergenceMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GenericMap.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmOutput.h"

// Datamodel for INPUT products
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsPeakCountConvergence.h"
// Datamodel for OUTPUT products
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/out/euc-test-le3-wl-twodmass-PeakCatalog.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

using boost::program_options::options_description;
using boost::program_options::variable_value;

using std::string;
using std::chrono::duration;
using std::chrono::high_resolution_clock;
using LE3_2D_MASS_WL_PEAK_COUNT::PeakParam;
using LE3_2D_MASS_WL_CARTESIAN::GenericMap;
using LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap;
using namespace LE3_2D_MASS_WL_UTILITIES;

// Input namespace and classes
using namespace dpd::le3::wl::twodmass::inp::paramspeakcountconvergence;

// Output namespace and classes
using namespace dpd::le3::wl::twodmass::out::peakcatalog;
//(class dpdTwoDMassPeakCatalog)

static Elements::Logging logger = Elements::Logging::getLogger(
        "PeakCountConvergence");

class PeakCountConvergence: public Elements::Program
{

public:

    options_description defineSpecificProgramOptions() override
    {

        options_description options { };
        // input parameter file
        options.add_options()("paramPeakConvergence",
                po::value<string>()->default_value(""),
                "input configuration file");
        // input working directory
        options.add_options()("workdir", po::value<string>(),
                "Work Directory Path");
        // input convergence map (cartesian)
        options.add_options()("inputConvMap",
                po::value<string>()->default_value(""),
                "input Convergence Map");
        // output file for Peak Catalog
        options.add_options()("outputPConvCatalog",
                po::value<string>()->default_value(""),
                "Output Peak Count Catalog in fits");
        options.add_options()("outputPConvCatalogXML",
                po::value<string>()->default_value("outputPConvCatalogXML.xml"),
                "Output Peak Count Catalog in XML");
        return options;
    }

    Elements::ExitCode mainMethod(std::map<std::string, variable_value>& args)
    override
    {
        /*
        ////////////////////////////////////////////////////////////////////////
        // start time & intro
        ////////////////////////////////////////////////////////////////////////
        auto start_ConvPeakTime = std::chrono::system_clock::now();
        logger.info("#");
        logger.info("# Entering mainMethod()");
        logger.info("#");

        ////////////////////////////////////////////////////////////////////////
        // get out time running ID for filename part
        ////////////////////////////////////////////////////////////////////////
        logger.info() << "Using " << omp_get_max_threads() << " thread(s) max.";

        ////////////////////////////////////////////////////////////////////////
        // Get the workdir and the data_dir
        ////////////////////////////////////////////////////////////////////////
        fs::path workdir { args["workdir"].as<string>() };
        fs::path datadir { workdir / "data" };
        fs::path fileParamPeakConv
                             { args["paramPeakConvergence"].as<std::string>() };
        fs::path inputConvMap { args["inputConvMap"].as<std::string>() };
        fs::path outputPConvCatalog
                               { args["outputPConvCatalog"].as<std::string>() };

        ////////////////////////////////////////////////////////////////////////
        // check parameter file exists
        ////////////////////////////////////////////////////////////////////////
        if (false == fs::is_regular(workdir / fileParamPeakConv))
        {
            logger.info() << "Parameter file not found or not in XML format...";
            return Elements::ExitCode::NOINPUT;
        }

        ////////////////////////////////////////////////////////////////////////
        // check whether parameter file (XML format) is for wavelet based peak
        // count Algorithm or not
        ////////////////////////////////////////////////////////////////////////
        PeakParam params;
        readPeakParamFile((workdir / fileParamPeakConv), params);

        ////////////////////////////////////////////////////////////////////////
        // Output filename
        ////////////////////////////////////////////////////////////////////////
        if ((outputPConvCatalog.string()).empty() == true)
        {
            outputPConvCatalog = fs::path(
                    "EUC_LE3_WL_PeakCatalog_Wavelet_" + getDateTimeString()
                            + ".fits");
        }

        ////////////////////////////////////////////////////////////////////////
        //Peak Count Wavelet algorithm
        ////////////////////////////////////////////////////////////////////////
        ConvergenceMap conv(workdir / inputConvMap);
        logger.info() << "PixelSize: " << conv.getPixelSize();

        WaveletPeakCount myPeakCounter(conv, params);
        myPeakCounter.savePeakCatalog((datadir / outputPConvCatalog).string());

        ////////////////////////////////////////////////////////////////////////
        // Generate the output XML product
        ////////////////////////////////////////////////////////////////////////
        auto out_xml_peak_convergence_catalog =
                args["outputPConvCatalogXML"].as<std::string>();
        // Create a DMOutput object
        DmOutput dm;
        std::string catfile = (workdir
                / fs::path(out_xml_peak_convergence_catalog)).string();
        fs::path outfile(catfile);
        fs::path catOut((workdir / outputPConvCatalog).native());
        dm.createPeakOutputXml(outfile, catOut);

        logger.info() << "DM output products created in: " << catfile;

        logger.info("Done!");

        ////////////////////////////////////////////////////////////////////////
        //End of Peak Count Convergence Program
        ////////////////////////////////////////////////////////////////////////
        auto ConvPeak_end = std::chrono::system_clock::now();
        std::chrono::duration<double> ConvPeak_elapsed_seconds = ConvPeak_end
                - start_ConvPeakTime;
        logger.info() << "elapsed time: " << ConvPeak_elapsed_seconds.count()
                << "s";
        logger.info("#");
        logger.info("# Exiting mainMethod()");
        logger.info("#");
        */
        return Elements::ExitCode::OK;
    }
};

MAIN_FOR(PeakCountConvergence)
