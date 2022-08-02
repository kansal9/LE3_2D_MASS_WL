/**
 * @file src/program/LE3_2D_MASS_WL_PeakCountSphere.cpp
 * @date 06/18/21
 * @author vanshika kansal
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
#include "LE3_2D_MASS_WL_UTILITIES/DmOutput.h"
#include "LE3_2D_MASS_WL_PEAK_COUNT/SphericalPeakCount.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalIO.h"

// Datamodel for INPUT products
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsPeakCountConvergence.h"
// Datamodel for OUTPUT products
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/out/euc-test-le3-wl-twodmass-PeakCatalog.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace Euclid::WeakLensing::TwoDMass;
using namespace Euclid::WeakLensing::TwoDMass::Spherical;
using namespace LE3_2D_MASS_WL_UTILITIES;

using boost::program_options::options_description;
using boost::program_options::variable_value;

using std::string;
using std::chrono::duration;
using std::chrono::high_resolution_clock;
using LE3_2D_MASS_WL_PEAK_COUNT::PeakParam;
using LE3_2D_MASS_WL_PEAK_COUNT::SphericalPeakCount;
using Euclid::WeakLensing::TwoDMass::Spherical::SphericalIO;

// Input namespace and classes
using namespace dpd::le3::wl::twodmass::inp::paramspeakcountconvergence;

// Output namespace and classes
using namespace dpd::le3::wl::twodmass::out::peakcatalog;
//(class dpdTwoDMassPeakCatalog)

static Elements::Logging logger = Elements::Logging::getLogger(
        "PeakCountSphere");

class PeakCountSphere: public Elements::Program
{

public:

    OptionsDescription defineSpecificProgramOptions() override
    {

        OptionsDescription options{ };

        // input working directory
        options.add_options()("workdir", po::value<string>(),
                "Work Directory Path");
        // input parameter file
        options.add_options()("paramPeakConvergence",
                po::value<string>()->default_value(""),
                "input configuration file");
        // input convergence healpix map
        options.add_options()("inputConvMap",
                po::value<string>()->default_value(""),
                "input Convergence Map in Healpix format");
        // output file for Peak Catalog (fits and xml)
        options.add_options()("outputPeakCatalog",
                po::value<string>()->default_value(""),
                "Output Peak Count Catalog in fits");
        options.add_options()("outputPeakCatalogXML",
                po::value<string>()->default_value("outputPConvCatalogXML.xml"),
                "Output Peak Count Catalog in XML");
        return options;
    }

    Elements::ExitCode mainMethod(std::map<std::string, VariableValue>& args)
    override
    {
        /*
        ////////////////////////////////////////////////////////////////////////
        // start time & intro
        ////////////////////////////////////////////////////////////////////////
        auto start_PeakTime = std::chrono::system_clock::now();

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
        fs::path fileParam { args["paramPeakConvergence"].as<std::string>() };
        fs::path inputConvMap { args["inputConvMap"].as<std::string>() };
        fs::path outputPeakCatalog
        { args["outputPeakCatalog"].as<std::string>() };

        ////////////////////////////////////////////////////////////////////////
        // check parameter file exists
        ////////////////////////////////////////////////////////////////////////
        if (false == fs::is_regular(workdir / fileParam))
        {
            logger.info() << "Parameter file not found or not in XML format...";
            return Elements::ExitCode::NOINPUT;
        }

        ////////////////////////////////////////////////////////////////////////
        // check whether parameter file (XML format) is for wavelet based peak
        // count algorithm or not
        ////////////////////////////////////////////////////////////////////////
        PeakParam params;
        readPeakParamFile((workdir / fileParam), params);

        ////////////////////////////////////////////////////////////////////////
        // Output filename
        ////////////////////////////////////////////////////////////////////////
        if ((outputPeakCatalog.string()).empty() == true)
        {
            outputPeakCatalog = fs::path(
                    "EUC_LE3_WL_PeakCatalog_Sphere_" + getDateTimeString()
                            + ".fits");
        }

        ////////////////////////////////////////////////////////////////////////
        // Object to Read/Write Spherical Map
        ////////////////////////////////////////////////////////////////////////
        std::pair<Healpix_Map<double>, Healpix_Map<double> > mapPair;
        if ((false == inputConvMap.string().empty()))
        {
            mapPair = readHealpixMap((workdir / inputConvMap).native());
        }
        else
        {
            return Elements::ExitCode::NOINPUT;
        }

        Healpix_Map<double>& mapE = mapPair.first;

        ////////////////////////////////////////////////////////////////////////
        //Peak Count algorithm
        ////////////////////////////////////////////////////////////////////////
        SphericalPeakCount peakCounter(mapE, params);
        peakCounter.findPeaksAtTheta((datadir / outputPeakCatalog).native());

        ////////////////////////////////////////////////////////////////////////
        // Generate the output XML product
        ////////////////////////////////////////////////////////////////////////
        auto out_xml_peak_catalog =
                                 args["outputPeakCatalogXML"].as<std::string>();
        // Create a DMOutput object
        DmOutput dm;
        std::string catfile =
                            (workdir / fs::path(out_xml_peak_catalog)).string();
        fs::path outfile(catfile);
        fs::path catOut((workdir / outputPeakCatalog).native());
        dm.createPeakOutputXml(outfile, catOut);
        logger.info() << "DM output products created in: " << catfile;

        ////////////////////////////////////////////////////////////////////////
        //End of Peak Count Program
        ////////////////////////////////////////////////////////////////////////
        auto Peak_end = std::chrono::system_clock::now();
        std::chrono::duration<double> Peak_elapsed_seconds = Peak_end
                - start_PeakTime;
        logger.info() << "elapsed time: " << Peak_elapsed_seconds.count()<< "s";
        logger.info("#");
        logger.info("# Exiting mainMethod()");
        logger.info("#");
        */
        return Elements::ExitCode::OK;
    }
};

MAIN_FOR(PeakCountSphere)
