/**
 * @file src/program/LE3_2D_MASS_WL_PatchesToSphere.cpp
 * @date 02/21/20
 * @author vkansal
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

#include <boost/program_options.hpp>
#include "ElementsKernel/ProgramHeaders.h"

#include "LE3_2D_MASS_WL_CARTESIAN/ShearMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/MatrixProcess.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianAlgoKS.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/PatchesToSphereAlgo.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GenericMap.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalIO.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmOutput.h"

#include <algorithm>
#include <cmath>
#include <vector>

using boost::program_options::options_description;
using boost::program_options::variable_value;

using std::string;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

using namespace LE3_2D_MASS_WL_UTILITIES;
using namespace LE3_2D_MASS_WL_CARTESIAN;

/*
 //TODO write these maps (output of Patches to Sphere Pipeline)
 Convergence Map => path position array (contains centers of the patches for each path) => DONE
 galaxy count for the map => DONE
 Noisy Map => DONE
 Denoised map => DONE
 SNR Maps
 MonteCarlo maps (get json as input from program LE3_2D_MASS_WL_CartesianMCMaps.cpp) and get output.
 */

static Elements::Logging logger = Elements::Logging::getLogger(
        "PatchesToSphere");

class PatchesToSphere: public Elements::Program
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

        // input file which contains list of Shear maps
        options.add_options()("inputShearMapList",
                po::value<string>()->default_value(""),
                "input Shear Maps list in json file format");

        // input parameter file
        options.add_options()("paramFile",
                po::value<string>()->default_value(""),
                "Input Parameter File in XML");

        // output file which contains list of centers for the Convergence maps
        options.add_options()("outputPatchesCenters",
                po::value<string>()->default_value(""),
                "output file containing the position of the centre of the "
                "patches used to compute the projections in fits format.");

        // output file in xml format for list of centers for the Convergence maps
        options.add_options()("outputPatchesCentersXML",
                po::value<string>()->default_value(""),
                "output file containing the position of the centre of the"
                "patches used to compute the projections in xml format.");

        // output file which contains Convergence maps values
        options.add_options()("outputHealpixConvergences",
                po::value<string>()->default_value(""),
                "Output list of Spherical convergence Maps (Noisy, Denoised and SNR) in json format.");

        // output file which contains Convergence maps values
        options.add_options()("outputHealpixConvergence",
                po::value<string>()->default_value(""),
                "Output convergence map E & B mode in BinTable (healpix) Format.");

        // output file in xml format for the Convergence maps
        options.add_options()("outputHealpixConvergenceXML",
                po::value<string>()->default_value(""),
                "Output XML product for patches to sphere convergence map.");

        options.add_options()("GalCountMap",
                po::value<string>()->default_value(""),
                "output Galaxy count Map which conatins number of Galaxies per"
                "pixel for each redshift bin in fits format");

        options.add_options()("GalCountMapXML",
                po::value<string>()->default_value(""),
                "output Galaxy count Map which conatins number of Galaxies per"
                "pixel for each redshift bin in XML format");

        return options;
    }

    Elements::ExitCode mainMethod(std::map<std::string, variable_value>& args)
    override
    {

        logger.info("#");
        logger.info("# Entering mainMethod()");
        logger.info("#");

        ////////////////////////////////////////////////////////////////////////
        // Get the workdir
        ////////////////////////////////////////////////////////////////////////
        fs::path workdir { args["workdir"].as<string>() };
        fs::path datadir { workdir / "data" };
        fs::path logPath(args["logdir"].as<string>());
        fs::path ParamFile { args["paramFile"].as<std::string>() };
        fs::path inputShearMaps { args["inputShearMapList"].as<std::string>() };
        fs::path outputPatchesCenters
        { args["outputPatchesCenters"].as<std::string>() };
        fs::path outputPatchesCentersXML
        { args["outputPatchesCentersXML"].as<std::string>() };
        fs::path outputHealpixConvergence
        { args["outputHealpixConvergence"].as<std::string>() };
        fs::path outputHealpixConvergenceXML
        { args["outputHealpixConvergenceXML"].as<std::string>() };
        fs::path GalCountMap { args["GalCountMap"].as<std::string>() };
        fs::path GalCountMapXML { args["GalCountMapXML"].as<std::string>() };
        fs::path outConvergenceMaps { };
        fs::path outputHealpixConvMaps
        { args["outputHealpixConvergences"].as<std::string>() };

        ////////////////////////////////////////////////////////////////////////
        // data to store map position (Ra & Dec) and Pixel values
        ////////////////////////////////////////////////////////////////////////
        std::vector<std::vector<double>> data;

        ////////////////////////////////////////////////////////////////////////
        // check parameter file exists
        ////////////////////////////////////////////////////////////////////////
        if (false == fs::is_regular(workdir / ParamFile))
        {
            logger.info()
                    << "Parameter file not found or not in XML format ....";
            return Elements::ExitCode::NOINPUT;
        }

        ////////////////////////////////////////////////////////////////////////
        if (outputPatchesCenters.string().empty())
        {
            outputPatchesCenters = fs::path(
                    "outputPatchesCenters_" + getDateTimeString() + ".fits");
        }

        ////////////////////////////////////////////////////////////////////////
        // set output filename (json)
        ////////////////////////////////////////////////////////////////////////
        if ((outputHealpixConvMaps.string()).empty())
        {
            outputHealpixConvMaps = fs::path(
                    "SphericalConvergenceMapsList_" + getDateTimeString()
                            + ".json");
        }

        ////////////////////////////////////////////////////////////////////////
        // check whether parameter file is of XML format or not and read
        ////////////////////////////////////////////////////////////////////////
        CartesianParam params;
        CatalogData dummy;
        readParameterFile((workdir / ParamFile), params, dummy);

        ////////////////////////////////////////////////////////////////////////
        // Object to perform Cartesian KS algorithm
        ////////////////////////////////////////////////////////////////////////
        CartesianAlgoKS CartesainAlgo(params);

        ////////////////////////////////////////////////////////////////////////
        // Object to use functions for patches to sphere algorithm
        ////////////////////////////////////////////////////////////////////////
        PatchesToSphereAlgo p2s(params);

        ////////////////////////////////////////////////////////////////////////
        // get NSide for output maps
        ////////////////////////////////////////////////////////////////////////
        int nside = params.getNside();

        ////////////////////////////////////////////////////////////////////////
        outConvergenceMaps = fs::path(
                "NoisySphericalConvergenceMapsList_" + getDateTimeString()
                        + ".json");

        ////////////////////////////////////////////////////////////////////////
        // get Noisy Spherical convergence map
        ////////////////////////////////////////////////////////////////////////
        //TODO In case patches needs to be saved then don't use this fuction
        // instead use MassMapping program in pipeline
        // Delete getNoisyConvergenceMaps or getDeNoisyConvergenceMaps functions
        // (even from the class dedicated to patches2sphere)

        p2s.getNoisyConvergenceMaps(workdir, inputShearMaps,
                outConvergenceMaps);

        p2s.getPatchData(outConvergenceMaps, workdir,
                (datadir / outputPatchesCenters).string(), data);

        Healpix_Map<double> mapE, mapB, mapGalCount;
        p2s.getHealpixFormatMap(data, mapE, mapB, mapGalCount);

        // set FITS filenames
        if ((outputHealpixConvergence.string()).empty())
        {
            outputHealpixConvergence = fs::path(
                    "EUC_LE3_WL_PatchesToSphere_NoisyConvergenceMap_NSide"
                            + std::to_string(nside) + "_" + getDateTimeString()
                            + ".fits");
        }
        if ((GalCountMap.string()).empty())
        {
            GalCountMap = fs::path(
                    "EUC_LE3_WL_PatchesToSphere_GalCount_NSide"
                            + std::to_string(nside) + "_" + getDateTimeString()
                            + ".fits");
        }

        data.clear();

        // Writing Fits files
        logger.info("# Writing noised healpix Maps");
        p2s.write_Map((datadir / outputHealpixConvergence).native(), mapE,
                "KAPPA_E");
        p2s.write_Map((datadir / outputHealpixConvergence).native(), mapB,
                "KAPPA_B");
        p2s.write_Map((datadir / GalCountMap).native(), mapGalCount,
                "GALCOUNT");

        ////////////////////////////////////////////////////////////////////////
        // get Noisy Spherical convergence map
        ////////////////////////////////////////////////////////////////////////
        if (params.getGaussStd() != 0)
        {
            outConvergenceMaps = fs::path(
                    "DenoisedSphericalConvergenceMapsList_"
                            + getDateTimeString() + ".json");
            p2s.getDeNoisyConvergenceMaps(workdir, inputShearMaps,
                    outConvergenceMaps);
            p2s.getPatchData(outConvergenceMaps, workdir,
                    (datadir / outputPatchesCenters).string(), data);

            Healpix_Map<double> DmapE, DmapB, DmapGalCount;
            p2s.getHealpixFormatMap(data, DmapE, DmapB, DmapGalCount);

            outputHealpixConvergence = fs::path(
                    "EUC_LE3_WL_PatchesToSphere_DeNoisedConvergenceMap_NSide"
                            + std::to_string(nside) + "_" + getDateTimeString()
                            + ".fits");
            GalCountMap = fs::path(
                    "EUC_LE3_WL_PatchesToSphere_GalCount_NSide"
                            + std::to_string(nside) + "_" + getDateTimeString()
                            + ".fits");

            data.clear();

            // Writing Fits files
            logger.info("# Writing Denoised healpix Maps");
            p2s.write_Map((datadir / outputHealpixConvergence).native(), DmapE,
                    "KAPPA_E");
            p2s.write_Map((datadir / outputHealpixConvergence).native(), DmapB,
                    "KAPPA_B");
            p2s.write_Map((datadir / GalCountMap).native(), DmapGalCount,
                    "GALCOUNT");
        }

        //TODO for NReSamples

        logger.info("#");
        logger.info("# Exiting mainMethod()");
        logger.info("#");

        return Elements::ExitCode::OK;
    }

};

MAIN_FOR(PatchesToSphere)
