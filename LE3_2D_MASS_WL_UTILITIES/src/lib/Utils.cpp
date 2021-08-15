/**
 * @file src/lib/Utils.cpp
 * @date 10/13/20
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

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include <deque>
#include <iostream>

using ST_DM_Schema::getDmSchemaFilePath;
static Elements::Logging logger = Elements::Logging::getLogger("Utils");

namespace Euclid {
 namespace WeakLensing {
  namespace TwoDMass {

 bool checkFileType(fs::path fileName, const std::vector<char>& sign){
  // open file in binary mode
  std::fstream file(fileName.string(), std::fstream::in | std::fstream::binary);
  bool result = true;
  // compare result against magic numbers
  for (size_t i=0; i<sign.size(); ++i) {
   char c;
   file.get(c);
   result &= (c == sign[i]);
  }
  file.close(); //close the file
  return result;
 }

 // get DateTime string
 std::string getDateTimeString() {
  // Get current time from the clock, using microseconds resolution
  const boost::posix_time::ptime now = 
        boost::posix_time::microsec_clock::local_time();

  // Get the time offset in current day
  const boost::posix_time::time_duration td = now.time_of_day();
  const long month= now.date().month();
  const long day= now.date().day();
  const long year= now.date().year();
  //
  // Extract hours, minutes, seconds and milliseconds.
  //
  // Since there is no direct accessor ".milliseconds()",
  // milliseconds are computed _by difference_ between total milliseconds
  // (for which there is an accessor), and the hours/minutes/seconds
  // values previously fetched.
  //
  const long hours        = td.hours();
  const long minutes      = td.minutes();
  const long seconds      = td.seconds();
  const long milliseconds = td.total_milliseconds() -
                              ((hours * 3600 + minutes * 60 + seconds) * 1000);
  char buf[40];
  sprintf(buf, "%ld%02ld%02ldT%02ld%02ld%02ld.%03ldZ", year, month, day,
        hours, minutes, seconds, milliseconds);

  return std::string(buf);
 }

// parse file to search for matching string
 bool fileHasField(fs::path fileName, const std::string& key){
  bool found  = false;
  std::fstream inFile((fileName.string()).c_str(), std::ifstream::in);
  // loop file
  while (true == inFile.good()) {
   std::string line;
   getline(inFile, line); // get line from file
   if (std::string::npos != line.find(key)) {
    found = true;
    break;
   }
  }
 return found;
 }

//split string (delimiter = ".") {example: le3.wl.2dmass.output.patchconvergence}
// example: le3.wl.2dmass.output.patchconvergence to "le3" "wl" "2dmass" "output" "patchconvergence"
std::vector<std::string> split_fits_id(const std::string& text) {
  std::vector<std::string> results;
  boost::split(results, text, boost::is_any_of("."));
/*
  boost::split(results, text, [](char c){
    return c == '.';
  });*/
  results.shrink_to_fit();
  return results;
}

std::vector<fs::path> read_filenames(fs::path workdir, fs::path inputJsonFile) {
	std::ifstream file_stream(((workdir / inputJsonFile).string()).c_str());
	std::stringstream str_stream;
    str_stream << file_stream.rdbuf();
    std::string contents = str_stream.str();
	const std::string separators = "[, ]\t\n\r\"";
	boost::trim_if(contents, boost::is_any_of(separators));
	std::vector<fs::path> filenames;
	boost::split(filenames, contents, boost::is_any_of(separators), boost::token_compress_on);
    filenames.shrink_to_fit();
	return filenames;
}

bool vecMinMax(std::vector <double>& vector, double *vecMin, double *vecMax) {
 *vecMin=*vecMax=vector[0];
/* for (size_t i=1; i<vector.size(); i++) {
  if (vector[i]< *vecMin) {
   *vecMin=vector[i];
  }
  if (vector[i]> *vecMax) {
   *vecMax=vector[i];
  }
 }*/
auto const max_iter = std::max_element(vector.begin(), vector.end());
assert(max_iter != vector.end());
*vecMax = *max_iter;
auto const min_iter = std::min_element(vector.begin(), vector.end());
assert(min_iter != vector.end());
*vecMin = *min_iter;
return true;
}
/*
bool checkFileIsEuclidize(const std::string& filename, std::vector<std::string>& colname){

 MefFile b(filename, MefFile::Permission::Read);

 std::vector<std::string> ext =b.access<BintableHdu>(1).readColumnNames();

 std::string opStr = getColName(colname, "CORRECTION");
 for(size_t i = 0; i < colname.size(); i++) {
   if (std::find(ext.begin(), ext.end(), colname[i]) != ext.end()) {
     logger.info() << "column name -- " << colname[i] << " --found";
   } else {
     logger.info() << "column name -- " << colname[i] << " --not found";
     if (!colname[i].compare(opStr)) {
      return true;
     } else {
       logger.info() << "Input catalogue is not Euclidize catalogue..";
       return false;
     }
   }
 }
return true;
}*/

std::vector<std::string> getcolumnNames (const std::string& catType) {
 const std::string g1 = "SHE_" + catType + "_G1";
 const std::string g2 = "SHE_" + catType + "_G2";
 const std::string w = "SHE_" + catType + "_WEIGHT";
 const std::string ra = "SHE_" + catType + "_UPDATED_RA";
 const std::string dec = "SHE_" + catType + "_UPDATED_DEC";
 const std::string phzCorrection = "PHZ_" + catType + "_CORRECTION";
 const std::string z = "PHZ_MEDIAN";

 std::vector<std::string> colname;
 colname.push_back(ra);
 colname.push_back(dec);
 colname.push_back(g1);
 colname.push_back(g2);
 colname.push_back(w);
 colname.push_back(z);
 colname.push_back(phzCorrection);

 return colname;
}

int getIndexCol(std::vector<std::string>& colname, std::string string) {
   int ind = -1;
   auto it = std::find(colname.begin(), colname.end(), string);
   if (it == colname.end()) {
    logger.info() << "not in colnames";
   } else {
    ind = std::distance(colname.begin(), it);
    logger.info() << "column name " << string << " index is: " << ind;
  }
  return ind;
}

std::string getColName(std::vector<std::string>& colname, std::string substring1, std::string substring2) {
  std::string name;
  for (size_t i = 0; i < colname.size(); i++) {
    if (colname[i].find(substring1) != std::string::npos || colname[i].find(substring2) != std::string::npos ) {
      name = colname[i];
    }
  }
  return name;
}

std::string getColName(std::vector<std::string>& colname, std::string substring) {
  std::string name;
  for (size_t i = 0; i < colname.size(); i++) {
    if (colname[i].find(substring) != std::string::npos) {
      name = colname[i];
    }
  }
  return name;
}

std::pair<Healpix_Map<double>, Healpix_Map<double> > readHealpixMap(const std::string &Map) {
  arr<double> myarrE;
  arr<double> myarrB;

  MefFile fitsFile(Map, MefFile::Permission::Read);
  const auto &ext = fitsFile.access<BintableHdu>(1);
  const auto records = ext.parseAllRecords<boost::any>();
  const auto Ttype1 = records.as<std::string>("TTYPE1");
  const auto Ttype2 = records.as<std::string>("TTYPE2");
  //read ordering from header
  std::string ordering = records.as<std::string>("ORDERING");
  //logger.info() << "ordering of map: " << std::string(ordering);
  if (((Ttype1.value).compare("GAMMA1") == 0) || ((Ttype1.value).compare("KAPPA_E") == 0) ||
       ((Ttype1.value).compare("PIXEL") == 0)) {
   const auto mE = fitsFile.access<BintableHdu>(1).readColumn<double>(Ttype1);
   myarrE.copyFrom(mE.vector());
  }
  if (((Ttype2.value).compare("GAMMA2") == 0) || ((Ttype2.value).compare("KAPPA_B") == 0) ||
       ((Ttype2.value).compare("WEIGHT") == 0)) {
   const auto mB = fitsFile.access<BintableHdu>(1).readColumn<double>(Ttype2);
   myarrB.copyFrom(mB.vector());
  }

  Healpix_Map<double> mapE, mapB;
  //mapE.Set(myarrE, RING);
  //mapB.Set(myarrB, RING);
  mapE.Set(myarrE, ordering=="RING" ? RING : NEST);
  mapB.Set(myarrB, ordering=="RING" ? RING : NEST);

  return std::pair<Healpix_Map<double>, Healpix_Map<double> > (mapE, mapB);
}

std::pair<double, double> getIndex2RaDec(Healpix_Map<double>& map, int index) {
 pointing ptg;
 ptg = map.pix2ang(index);
 //logger.info() << "dec: " << -( ptg.theta - M_PI*double(0.5) )/deg2rad;
 //logger.info() << "ra: " << (ptg.phi/deg2rad + rotPhi);
 double dec = (-( ptg.theta - M_PI*double(0.5) )/deg2rad);
 double ra = (ptg.phi/deg2rad + rotPhi);
 return std::pair<double, double> (ra, dec);
}

} /* namespace TwoDMass */
} /* namespace WeakLensing */
} /* namespace Euclid */
