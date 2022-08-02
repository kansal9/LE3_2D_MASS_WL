/**
 * @file src/program/LE3_2D_MASS_WL_PeakCountShear.cpp
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

#include "ST_DataModelBindings/dpd/le3/wl/twodmass/out/euc-test-le3-wl-twodmass-PeakCatalog.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsPeakCountMassAperture.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergencePatch.h"

#include <boost/program_options.hpp>

#include "ElementsKernel/ProgramHeaders.h"

#include "LE3_2D_MASS_WL_UTILITIES/PatchDef.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GenericMap.h"
#include "LE3_2D_MASS_WL_PEAK_COUNT/PeakParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianAlgoKS.h"
#include "LE3_2D_MASS_WL_UTILITIES/NoisyCatalogData.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmOutput.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

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
using LE3_2D_MASS_WL_CARTESIAN::CartesianParam;

using namespace LE3_2D_MASS_WL_CARTESIAN;
using namespace LE3_2D_MASS_WL_PEAK_COUNT;
using namespace LE3_2D_MASS_WL_UTILITIES;

using namespace dpd::le3::wl::twodmass::inp::paramspeakcountmassap;
using namespace dpd::le3::wl::twodmass::inp::paramsconvergencepatch;
using namespace dpd::le3::wl::twodmass::out::peakcatalog;

static Elements::Logging logger = Elements::Logging::getLogger(
        "PeakCountShear");

class PeakCountShear: public Elements::Program
{

public:

    options_description defineSpecificProgramOptions() override
    {

        options_description options
        { };
        // input parameter file for peak count
        options.add_options()("paramPeakMassAperture",
                po::value<string>()->default_value(""),
                "input configuration/Parameter file");
        // input parameter file to get shear map
        options.add_options()("paramPatchMap",
                po::value<string>()->default_value(""),
                "input configuration/Param file to extract shear map");
        // input working directory
        options.add_options()("workdir", po::value<string>(),
                "Work Directory Path");
        //input catalog
        options.add_options()("inputShearCatalog",
                po::value<string>()->default_value(""), "input Shear catalog");
        // output Peak Catalog
        options.add_options()("outputPeakMACatalog",
                po::value<string>()->default_value(""),
                "Output Peak Count Catalog in fits");
        options.add_options()("outputPeakMACatalogXML",
                po::value<string>()->default_value(
                        "outputPeakMACatalogXML.xml"),
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
        auto start_MAPeakTime = std::chrono::system_clock::now();

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
        fs::path fileParamPeakMA
        { args["paramPeakMassAperture"].as<std::string>() };
        fs::path logPath(args["logdir"].as<string>());
        fs::path inputShearCat { args["inputShearCatalog"].as<std::string>() };
        fs::path filePatchParam { args["paramPatchMap"].as<std::string>() };
        fs::path outputPeakMACatalog
        { args["outputPeakMACatalog"].as<std::string>() };

        ////////////////////////////////////////////////////////////////////////
        // check parameter file exists
        ////////////////////////////////////////////////////////////////////////
        if ((false == fs::is_regular(workdir / fileParamPeakMA))
         || (false == fs::is_regular(workdir / filePatchParam)))
        {
            logger.info() << "Parameter file not found or not in XML format...";
            return Elements::ExitCode::NOINPUT;
        }

        ////////////////////////////////////////////////////////////////////////
        // check whether parameter file (XML format) is for Mass Aperture Peak
        // Count Algorithm or not
        ////////////////////////////////////////////////////////////////////////
        PeakParam params;
        readPeakParamFile((workdir / fileParamPeakMA), params);

        ////////////////////////////////////////////////////////////////////////
        // Reading input shear catalog
        ////////////////////////////////////////////////////////////////////////
        CatalogData cat;
        logger.info("Reading Input Shear Catalog");
        cat.getCatalogData(workdir, inputShearCat);

        ////////////////////////////////////////////////////////////////////////
        // Variables need in case of Patches to sphere Parameter file
        ////////////////////////////////////////////////////////////////////////
        // TODO ?

        ////////////////////////////////////////////////////////////////////////
        // check whether xml parameter file is for Patch based map or not
        ////////////////////////////////////////////////////////////////////////
        CartesianParam cartesianParam;
        readParameterFile((workdir / filePatchParam), cartesianParam, cat);
        //TODO read patches to sphere parameter file and get the peaks for all patches
        //readParameterFile ((workdir/filePatchParam), Patchparams, catRamin, catRamax, catDecmin, catDecmax);

        ////////////////////////////////////////////////////////////////////////
        // Object to fetch Shear Map
        ////////////////////////////////////////////////////////////////////////
        CartesianAlgoKS CartesainAlgo(cartesianParam);

        double zMin, zMax;
        vecMinMax(cat["z"], zMin, zMax);

        ////////////////////////////////////////////////////////////////////////
        // Creating Patch of shear Map from Input shear catalog
        ////////////////////////////////////////////////////////////////////////
        for (size_t it = 0; it < cartesianParam.getNPatches(); it++)
        {
            auto& patch = cartesianParam.getPatches()[it];
            logger.info("creating shear Map");
            double ramin, decmin;
            ramin = patch.getRaMin();
            decmin = patch.getDecMin();
            //PatchDef m_CB(ramin, cartesianParam.getRaMax(ramin), decmin,
            //        cartesianParam.getDecMax(decmin), zMin, zMax);
            fs::path shearMap = fs::path(
                    "EUC_LE3_WL_ShearMap_" + getDateTimeString() + ".fits");
            CartesainAlgo.extractShearMap(datadir / shearMap, cat, patch);

            ////////////////////////////////////////////////////////////////////
            // Randomising shear catalog
            ////////////////////////////////////////////////////////////////////
            NoisyCatalogData random;
            CatalogData noisyCat = random.getNoisyCatalog(cat);

            ////////////////////////////////////////////////////////////////////
            // Creating Patch of shear Map from randomised shear
            /////////////////////////////////////////////////////////////////////
            logger.info("creating Noisy shear Map");

            fs::path NoisyshearMap = fs::path(
                    "EUC_LE3_WL_NoisyShearMap_" + getDateTimeString()
                            + ".fits");
            CartesainAlgo.extractShearMap(datadir/NoisyshearMap,noisyCat, patch);

            ////////////////////////////////////////////////////////////////////
            // Output filename
            ////////////////////////////////////////////////////////////////////
            if ((outputPeakMACatalog.string()).empty())
            {
                outputPeakMACatalog = fs::path(
                        "EUC_LE3_WL_PeakCatalogMA_" + getDateTimeString()
                                + ".fits");
            }

            ////////////////////////////////////////////////////////////////////
            //Peak Count Mass Aperture algorithm
            ////////////////////////////////////////////////////////////////////
            ShearMap shear(datadir / shearMap);

            ShearMap noisedShear(datadir / NoisyshearMap);

            MassAperturePeakCount myPeakCounter(shear, noisedShear, params);
            myPeakCounter.saveMAPeakCatalog(
                                      (datadir / outputPeakMACatalog).string());

            //TODO write single peak catalogue (canbe done after adding dependency to EL_FtsIO new version)
        }

        ////////////////////////////////////////////////////////////////////////
        // Generate the output XML product
        ////////////////////////////////////////////////////////////////////////
        auto out_xml_peak_ma_catalog =
                args["outputPeakMACatalogXML"].as<std::string>();
        // Create a DMOutput object
        DmOutput dm;
        std::string catfile =
                (workdir / fs::path(out_xml_peak_ma_catalog)).string();
        fs::path outfile(catfile);
        fs::path catOut((workdir / outputPeakMACatalog).native());
        dm.createPeakOutputXml(outfile, catOut);
        logger.info() << "DM output products created in: " << catfile;

        ////////////////////////////////////////////////////////////////////////
        //End of Peak Count MA Program
        ////////////////////////////////////////////////////////////////////////
        auto MAPeak_end = std::chrono::system_clock::now();
        std::chrono::duration<double> MAPeak_elapsed_seconds = MAPeak_end
                - start_MAPeakTime;
        logger.info() << "elapsed time: "<<MAPeak_elapsed_seconds.count()<<"s";
        logger.info("#");
        logger.info("# Exiting mainMethod()");
        logger.info("#");
        */
        return Elements::ExitCode::OK;
    }
};

MAIN_FOR(PeakCountShear)
