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
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>
#include <algorithm>
#include <cstdlib>

// Datamodel for INPUT products
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergencePatch.h"

#include <boost/program_options.hpp>
#include "ElementsKernel/ProgramHeaders.h"

#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"
#include "LE3_2D_MASS_WL_CARTESIAN/SplitCatalog.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_UTILITIES/ReadCatalog.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

using boost::program_options::options_description;
using boost::program_options::variable_value;
using std::chrono::duration;
using std::chrono::high_resolution_clock;
using std::string;
using std::vector;

using namespace Euclid::WeakLensing::TwoDMass;
using LE3_2D_MASS_WL_CARTESIAN::CartesianParam;
using LE3_2D_MASS_WL_SPHERICAL::SphericalParam;

static Elements::Logging logger = Elements::Logging::getLogger("LE3_2D_MASS_WL_CatalogSplitter");

/**
 *  @brief    Template to sort indices of a vector
 *  @param    v <vector<T> >, vector to sort
 *  @return   vector of sorted indices
*/
template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) {
    idx[i] = i;
  }
  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {
   return v[i1] < v[i2];
  });
  return idx;
}

void getSortedData(std::vector<std::vector<double> >& Data) {
   std::vector <double> sorted_ra, sorted_dec, sorted_kappa, sorted_g1, sorted_g2, sorted_z, sorted_weight;
   for (auto i : sort_indexes(Data[5])) {
    sorted_ra.push_back(Data[0][i]);
    sorted_dec.push_back(Data[1][i]);
    sorted_z.push_back(Data[5][i]);
    sorted_kappa.push_back(Data[2][i]);
    sorted_g1.push_back(Data[3][i]);
    sorted_g2.push_back(Data[4][i]);
    sorted_weight.push_back(Data[6][i]);
   }
   Data.clear();
   Data.push_back(sorted_ra);
   Data.push_back(sorted_dec);
   Data.push_back(sorted_kappa);
   Data.push_back(sorted_g1);
   Data.push_back(sorted_g2);
   Data.push_back(sorted_z);
   Data.push_back(sorted_weight);
}

// This makes the sort be according to redshift column and ascending
/*bool sortFunc( const vector<double>& z1,
           const vector<double>& z2 ) {
 return z1[5] < z2[5];
 }*/

std::string readParamFile (const boost::filesystem::path& parameterFile,
     LE3_2D_MASS_WL_CARTESIAN::CartesianParam &params, LE3_2D_MASS_WL_SPHERICAL::SphericalParam &Sparams,
     int& nbZbin, long& balancedBins) {
  double placeHolder = 0.0;
  std::string m_split;
   // Case 1: check parameter file is of XML format
   if (true == checkFileType(parameterFile.native(), Euclid::WeakLensing::TwoDMass::signXML)) {
    // read XML, get parameter file type
    if (true == fileHasField(parameterFile.native(), "DpdTwoDMassParamsConvergencePatch")) {
      logger.info()<<"Parameter file is for Convergence Patch..";
      params.ReadConvPatchXMLFile (parameterFile.native());
      nbZbin = params.getnbZBins();
      balancedBins = params.get_BalancedBins();
      m_split = "Cartesian";
    }
    if (true == fileHasField(parameterFile.native(), "DpdTwoDMassParamsConvergencePatchesToSphere")) {
      logger.info()<<"Parameter file is for Convergence Patches to Sphere..";
      params.readConvPatchesToSphereXMLFile (parameterFile.native(), placeHolder, placeHolder,
                                                placeHolder, placeHolder);
      nbZbin = params.getnbZBins();
      balancedBins = params.get_BalancedBins();
      m_split = "Cartesian";
    }
    if (true == fileHasField(parameterFile.native(), "DpdTwoDMassParamsConvergenceSphere")) {
      logger.info()<<"Parameter file is for Convergence for Sphere..";
      Sparams.getConvergenceSphereParam (parameterFile.native());
      nbZbin = Sparams.getNbins();
      balancedBins = Sparams.getBalancedBins();
      m_split = "Spherical";
   }
  }
 return m_split;
}

class LE3_2D_MASS_WL_CatalogSplitter : public Elements::Program {

public:

  options_description defineSpecificProgramOptions() override {
  
    options_description options {};
   // input working directory
   options.add_options()
   ("workdir", po::value<string>()->default_value(""), "Working Directory Path");

   // input log directory
   options.add_options()
   ("logdir", po::value<string>()->default_value(""), "log Directory Path");

   // input shear catalog file
   options.add_options()
   ("inputCatalog", po::value<string>()->default_value(""), "input shear Catalog");

   // input parameter file
   options.add_options()
   ("paramFile", po::value<string>()->default_value(""), "Input Parameter File to split input catalog");

   // output product (list of sub-catalogs)
   options.add_options()
//   ("Sub_Catalogs", po::value<std::vector<std::string> >()->default_value(std::vector<std::string>(), ""),
//                    "Sub-Catalog Name");
   ("Sub_Catalogs", po::value<string>()->default_value(""), "Sub-Catalog Name in txt/jason file");

    return options;
  }

  Elements::ExitCode mainMethod(std::map<std::string, variable_value>& args) override {
    //Start time and intro
    auto start_SplitterTime = std::chrono::system_clock::now();
    logger.info("#");
    logger.info("# Entering mainMethod()");
    logger.info("#");

////////////////////////////////////////////////////////////////////////////////////////////////////////
  // get out time running ID for filename part
////////////////////////////////////////////////////////////////////////////////////////////////////////
//  logger.info("Run ID: " + getDateTimeString());

  logger.info() << "Using " << omp_get_max_threads() << " thread(s) max.";

////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Get the workdir, the datadir and manage Input output
////////////////////////////////////////////////////////////////////////////////////////////////////////
   fs::path workdir {args["workdir"].as<string>()};
   fs::path datadir {workdir / "data"};
   fs::path Input_Catalog {args["inputCatalog"].as<string>()};
   fs::path parameterFile {args["paramFile"].as<std::string>()};
   fs::path logPath (args["logdir"].as<string>());
   std::vector<std::string> SubCatFilenames;// = args["Sub_Catalogs"].as<std::vector<std::string> >();
////////////////////////////////////////////////////////////////////////////////////////////////////////
  // check parameter file exists
////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (false == fs::is_regular(workdir /parameterFile)) {
   logger.info()<< "Parameter file not found or not in XML format ....";
   return Elements::ExitCode::NOINPUT;
  }
  int nbZbin;
  long m_balancedBins;
  std::string m_split;
////////////////////////////////////////////////////////////////////////////////////////////////////////
  // check whether parameter file is of XML format
////////////////////////////////////////////////////////////////////////////////////////////////////////
  CartesianParam params;
  SphericalParam Sparams;

  m_split = readParamFile ((workdir /parameterFile), params, Sparams, nbZbin, m_balancedBins);

 if (nbZbin <= 1) {
    SubCatFilenames.push_back(Input_Catalog.string());
 }

 if (nbZbin > 1) {
  //logger.info("Splitting input Catalog into sub-catalogs");
////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Variable to save input catalog name
////////////////////////////////////////////////////////////////////////////////////////////////////////
   std::string inputCatalog;
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Reading input catalog
////////////////////////////////////////////////////////////////////////////////////////////////////////
  ReadCatalog read;
  std::vector<std::vector<double> > catData;
  std::string m_catType="LENSMC";

  read.readShearCatalog(workdir, Input_Catalog, catData);
  m_catType = read.getCatalogReadType();

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// If redshift is 0 for all objects in catalog and balanced bin is false, then set balanced bin to true
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  double zMin, zMax;
  vecMinMax(catData[5], &zMin, &zMax);
  if (zMin == 0 && zMax == 0){
  //if (zMin == zMax){
   logger.info() << "ERROR: The redshift values for all objects inside catalog is 0.";
   return Elements::ExitCode::OK;
  }
  for (int i =0; i<nbZbin; i++) {
     std::string subCatalogName =  (fs::path ("EUC_LE3_WL_CatZ_subcatalogue_Z0" +
                                 std::to_string(i)+ "_" + getDateTimeString() + ".fits")).string();
     SubCatFilenames.push_back(subCatalogName);
  }
  if (m_balancedBins == 1) {
    getSortedData(catData);
  // Do the sorting according to redshift
   // std::sort(catData.begin(), catData.end(), sortFunc);
   //std::stable_sort(catData.begin(), catData.end(), sortFunc);
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Splitting Input catalog
////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (m_split == "Cartesian") {
    LE3_2D_MASS_WL_CARTESIAN::SplitCatalog split(catData, params, m_catType);
    split.writeSubCatalogs(datadir, SubCatFilenames);
  }
  if (m_split == "Spherical") {
    LE3_2D_MASS_WL_CARTESIAN::SplitCatalog split(catData, Sparams, m_catType);
    split.writeSubCatalogs(datadir, SubCatFilenames);
  }
catData.clear();
 }
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Generate the jason/txt product
////////////////////////////////////////////////////////////////////////////////////////////////////////
  fs::path output_filename = args["Sub_Catalogs"].as<string>();
////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (true == output_filename.string().empty()) {
      output_filename = fs::path("SubCatalogs_" + getDateTimeString() + ".json");
  }
  std::ofstream outfile ((workdir /output_filename).string());
  outfile << "[";
  size_t i = 0;
  while (i != SubCatFilenames.size()) {
    //fs::path subCatName(SubCatFilenames[i]);
    //outfile << subCatName.filename();
    outfile <<"\"" << SubCatFilenames[i] << "\"";
    i++;
    if (i<SubCatFilenames.size()) {
      outfile << ", ";
    }
  }
  outfile << "]";
 //}
////////////////////////////////////////////////////////////////////////////////////////////////////////
 //End of Splitter
////////////////////////////////////////////////////////////////////////////////////////////////////////
   auto end_SplitterTime = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = end_SplitterTime-start_SplitterTime;
   //std::time_t end_time = std::chrono::system_clock::to_time_t(end_SplitterTime);
   //logger.info() << "finished Program at " << std::ctime(&end_time);
   logger.info() << "elapsed time: " << elapsed_seconds.count() << "s";

    //logger.info("#");
    logger.info("# Exiting mainMethod()");
    logger.info("#");

    return Elements::ExitCode::OK;
  }

};

MAIN_FOR(LE3_2D_MASS_WL_CatalogSplitter)
