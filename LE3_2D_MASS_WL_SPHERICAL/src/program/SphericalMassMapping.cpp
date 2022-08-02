/**
 * @file src/program/LE3_2D_MASS_WL_SphericalMassMapping.cpp
 * @date 11/30/19
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

#include <boost/property_tree/ptree.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include "ElementsKernel/ProgramHeaders.h"

#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalIO.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalInpainting.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalUtils.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphMassMapping.h"

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmOutput.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"

// Datamodel for INPUT products
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergenceSphere.h"

// Datamodel for OUTPUT products
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/out/euc-test-le3-wl-twodmass-ConvergenceSphere.h"

using std::chrono::duration;
using std::chrono::high_resolution_clock;
using std::string;
using std::vector;

using LE3_2D_MASS_WL_SPHERICAL::SphMassMapping;
using LE3_2D_MASS_WL_SPHERICAL::SphericalParam;
using LE3_2D_MASS_WL_SPHERICAL::SphericalInpainting;
using Euclid::WeakLensing::TwoDMass::Spherical::SphericalIO;

// Input namespace and classes
using namespace dpd::le3::wl::twodmass::inp::paramsconvergencesphere;

// Output namespace and classes
using namespace dpd::le3::wl::twodmass::out::convergencesphere;

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using boost::program_options::options_description;
using boost::program_options::variable_value;
using namespace Euclid::WeakLensing::TwoDMass;
using namespace Euclid::WeakLensing::TwoDMass::Spherical;
using namespace LE3_2D_MASS_WL_UTILITIES;

static Elements::Logging logger = Elements::Logging::getLogger(
        "SphericalMassMapping");

/**
 *  @brief    function to perform Spherical mass mapping
 *  @param    mapPair, Gamma1 and Gamma2 in healpix format
 *  @param    kappaPair, kappaE and kappaB in healpix format (output)
 *  @param    SphericalParam, object that points to parameters
 *  @return   None
 */

void getSphericalConvergenceMap(
        std::pair<Healpix_Map<double>, Healpix_Map<double> >& mapPair,
        std::pair<Healpix_Map<double>, Healpix_Map<double> >& kappaPair,
        SphericalParam &SphericalParam)
{

    SphMassMapping mapping;

    kappaPair = mapping.create_SheartoConvMap(mapPair.first, mapPair.second);

    // Perform Inpainting

    if (SphericalParam.getNInpaint() > 0)
    {

        SphericalInpainting Inpainting(mapPair, kappaPair, SphericalParam);
        kappaPair = Inpainting.performInpainting();

    }

}

class SphericalMassMapping: public Elements::Program
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

        // input parameter file
        options.add_options()("paramFile",
                po::value<string>()->default_value(""),
                "Input Parameter File in XML");

        // input product (shear Map)
        options.add_options()("inGamma", po::value<string>()->default_value(""),
                "input shear Map E and B mode in fits format ");

        // input product (Density Map)
        options.add_options()("inGalCount",
                po::value<string>()->default_value(""),
                "Galaxy count Map which conatins number of Galaxies per pixel in fits format ");

        // input product (convergence Map)
        options.add_options()("inKappa", po::value<string>()->default_value(""),
                "input kappa Map E and B mode in fits format ");

        // output json file for convergance Maps
        options.add_options()("outputKappa",
                po::value<string>()->default_value("SphereConvergence.json"),
                "Output convergence map in json format that stores convergence E & B mode in BinTable Format");
        // output FITS file for Shear Map E and B mode
        options.add_options()("outputShear",
                po::value<string>()->default_value(""),
                "Output convergence map E mode in BinTable Format");
        // output XML for convergance Map E and B mode
        options.add_options()("outputKappaXML",
                po::value<string>()->default_value("outputKappaXML.xml"),
                "Output XML product for Spherical convergence map");
        return options;
    }

    Elements::ExitCode mainMethod(std::map<std::string, variable_value>& args)
    override
    {

        // start time & intro
        auto SphMassMapping_start = std::chrono::system_clock::now();

        //logger.info("#");
        logger.info("# Entering mainMethod()");
        // logger.info("#");
        ////////////////////////////////////////////////////////////////////////
        // get out time running ID for filename part
        ////////////////////////////////////////////////////////////////////////
        //logger.info("Run ID: " + getDateTimeString());
        logger.info() << "Using " << omp_get_max_threads() << " thread(s) max.";

        ////////////////////////////////////////////////////////////////////////
        // Get the workdir and manage Input output
        ////////////////////////////////////////////////////////////////////////
        fs::path workdir { args["workdir"].as<string>() };
        fs::path ParamFile { args["paramFile"].as<std::string>() };
        fs::path inputGamma { args["inGamma"].as<string>() };
        fs::path inputGalCount { args["inGalCount"].as<string>() };
        fs::path inputKappa { args["inKappa"].as<string>() };
        fs::path outputKappa { args["outputKappa"].as<string>() };
        fs::path outputGamma { args["outputShear"].as<string>() };

        fs::path datadir { workdir / "data" };
        fs::path outputKappaFits { };

        ////////////////////////////////////////////////////////////////////////
        // check parameter file and input file (Gamma or Kappa Map) exists
        ////////////////////////////////////////////////////////////////////////
        if (false == fs::is_regular(workdir / ParamFile))
        {
            logger.info()
                    << "Parameter file not found or not in XML format ....";
            return Elements::ExitCode::NOINPUT;
        }

        ////////////////////////////////////////////////////////////////////////
        // check whether parameter file is of XML format or .ini format
        ////////////////////////////////////////////////////////////////////////
        SphericalParam SphParam;
        CatalogData dummy;
        readSphericalParameterFile((workdir / ParamFile), SphParam, dummy);

        ////////////////////////////////////////////////////////////////////////
        // Object to Read/Write Spherical Map
        ////////////////////////////////////////////////////////////////////////
        SphericalIO SphericalIO(SphParam);
        std::pair<Healpix_Map<double>, Healpix_Map<double> > mapPair;
        if ((false == inputKappa.string().empty()))
        {
            mapPair = readHealpixMap((workdir / inputKappa).native());
        }
        else
        {
            mapPair = readHealpixMap((workdir / inputGamma).native());
        }

        Healpix_Map<double> mapE = mapPair.first;
        Healpix_Map<double> mapB = mapPair.second;

        ////////////////////////////////////////////////////////////////////////
        // Creating Shear Map from Convergence Map (inverse mass mapping KS)
        ////////////////////////////////////////////////////////////////////////
        if ((false == inputKappa.string().empty()))
        {
            logger.info("# Running Spherical inverse KS Mass Mapping");

            SphMassMapping mapping;
            std::pair<Healpix_Map<double>, Healpix_Map<double> > k2shearPair =
                    mapping.create_ConvtoShearMap(mapE, mapB);

            ////////////////////////////////////////////////////////////////////
            // Generate the Kappa to shear FITS filenames
            ////////////////////////////////////////////////////////////////////

            if (outputGamma.string().empty())
            {
                outputGamma = fs::path(
                        "EUC_LE3_WL_K2Shear_" + getDateTimeString() + ".fits");
            }
            ////////////////////////////////////////////////////////////////////
            // Writing Kappa 2 shear Fits files
            ////////////////////////////////////////////////////////////////////
            logger.info(
                    "# Writing Spherical Gamma Map got from Convergence Map");

            SphericalIO.write_Map((datadir / outputGamma).native(),
                    k2shearPair.first, "GAMMA1");
            SphericalIO.write_Map((datadir / outputGamma).native(),
                    k2shearPair.second, "GAMMA2");
        }

        ////////////////////////////////////////////////////////////////////////
        // Creating Convergence Map from Shear Map (Perform mass mapping KS)
        ////////////////////////////////////////////////////////////////////////
        logger.info("# Running Spherical Mass Mapping");

        std::pair<Healpix_Map<double>, Healpix_Map<double> > kappaPair;

        std::ofstream outfile;
        outfile.open((workdir / outputKappa).string(), std::ios_base::app);
        outfile << "[";

        ////////////////////////////////////////////////////////////////////////
        // Get Noisy Kappa Maps
        ////////////////////////////////////////////////////////////////////////
        getSphericalConvergenceMap(mapPair, kappaPair, SphParam);

        // Generate the output FITS filenames
        outputKappaFits = fs::path(
                "EUC_LE3_WL_NoisySphereConvergence_" + getDateTimeString()
                        + ".fits");

        // Writing output Fits files
        logger.info("# Writing Spherical Noisy Convergence Map");
        SphericalIO.write_Map((datadir / outputKappaFits).native(),
                kappaPair.first, "KAPPA_E");
        SphericalIO.write_Map((datadir / outputKappaFits).native(),
                kappaPair.second, "KAPPA_B");
        outfile << outputKappaFits.filename();

        ////////////////////////////////////////////////////////////////////////
        // Get De-Noisy Kappa Maps
        ////////////////////////////////////////////////////////////////////////

        if (fabs(SphParam.getGaussStd()) > 0.001)
        {
            applyGaussianFilter_hp(mapE, SphParam.getGaussStd()); //this will update contents by applying filter
            applyGaussianFilter_hp(mapB, SphParam.getGaussStd()); //this will update contents by applying filter
            outputKappaFits.clear();
            mapPair = std::make_pair(mapE, mapB);
            getSphericalConvergenceMap(mapPair, kappaPair, SphParam);
            // Generate the output FITS filenames
            outputKappaFits = fs::path(
                    "EUC_LE3_WL_DenoisedSphereConvergence_"
                            + getDateTimeString() + ".fits");
            outfile << ",";
            // Writing output Fits files
            logger.info("# Writing Spherical DeNoised Convergence Map");
            SphericalIO.write_Map((datadir / outputKappaFits).native(),
                    kappaPair.first, "KAPPA_E");
            SphericalIO.write_Map((datadir / outputKappaFits).native(),
                    kappaPair.second, "KAPPA_B");
            outfile << outputKappaFits.filename();
        }

        outfile << "]";
        outfile.close();

        ////////////////////////////////////////////////////////////////////////
        // Generate the output XML product
        ////////////////////////////////////////////////////////////////////////
        fs::path out_xml_product_file_name =
                args["outputKappaXML"].as<string>();
        SphericalIO.writeSphericalXMLfile(workdir, outputKappa, inputGalCount,
                out_xml_product_file_name);

        logger.info("Done!");

        ////////////////////////////////////////////////////////////////////////
        // End of the Program and Gives the final time of execution
        ////////////////////////////////////////////////////////////////////////
        auto SphMassMapping_end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = SphMassMapping_end
                - SphMassMapping_start;
        logger.info() << "elapsed time: " << elapsed_seconds.count() << "s";

        logger.info("# Exiting mainMethod()");
        logger.info("#");

        return Elements::ExitCode::OK;
    }

};

MAIN_FOR(SphericalMassMapping)
