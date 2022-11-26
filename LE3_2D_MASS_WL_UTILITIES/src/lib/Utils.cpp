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

//using ST_DM_Schema::getDmSchemaFilePath;
static Elements::Logging logger = Elements::Logging::getLogger("Utils");

namespace LE3_2D_MASS_WL_UTILITIES {

bool checkFileType(fs::path fileName, const std::vector<char>& sign)
{
    // open file in binary mode
    std::fstream file(fileName.string(),std::fstream::in|std::fstream::binary);
    // compare result with signature
    char c = (char) 0;
    bool result = true;
    for (size_t i = 0; i < sign.size(); ++i)
    {
        file.get(c);
        result &= (c == sign[i]);
    }
    file.close();
    return result;
}

// get DateTime string
std::string getDateTimeString()
{
    // Get current time from the clock, using microseconds resolution
    const pt::ptime now = pt::microsec_clock::local_time();
    return to_iso_string(now);
}

fs::path getFitsFilenameFromBase(const std::string& base)
{
    return fs::path(base + "_" + getDateTimeString() + ".fits");
}

// parse file to search for matching string
bool fileHasField(fs::path fileName, const std::string& key)
{
    bool found = false;
    std::fstream inFile((fileName.string()).c_str(), std::ifstream::in);
    std::string line;
    // loop file lines
    while (inFile.good())
    {
        getline(inFile, line);
        if (std::string::npos != line.find(key))
        {
            found = true;
            break;
        }
    }
    return found;
}

std::vector<std::string> splitString(const std::string& text,
                                     const std::string& split_string)
{
    std::vector<std::string> results;
    boost::split(results, text, boost::is_any_of(split_string));
    results.shrink_to_fit();
    return results;
}

std::vector<fs::path> readFilenamesInJson(fs::path inputJsonFile)
{
    std::ifstream file_stream(inputJsonFile.string().c_str());
    std::stringstream str_stream;
    str_stream << file_stream.rdbuf();
    std::string contents = str_stream.str();
    const std::string separators = "[, ]\t\n\r\"";
    boost::trim_if(contents, boost::is_any_of(separators));
    std::vector<fs::path> filenames;
    boost::split(filenames, contents, boost::is_any_of(separators),
            boost::token_compress_on);
    filenames.shrink_to_fit();
    return filenames;
}

std::vector<fs::path> readFilenames(fs::path inputFile)
{
    std::vector<fs::path> filenames;
    // Case: file does not exist: return empty vector
    if (!boost::filesystem::exists(inputFile))
    {
        return filenames;
    }
    // Case: file is in fits, fill vector with only one entry
    if (checkFileType(inputFile.native(), signFITS))
    {
        filenames.push_back(inputFile);
    }
    // Case: assume file is in json, open and get vector of filenames in json
    else
    {
        filenames = readFilenamesInJson(inputFile);
    }
    return filenames;
}

std::vector<std::pair<std::string, std::string>> getShearColNamesAndProxy(
        const std::string& shearType)
{
    std::vector<std::pair<std::string, std::string>> colnames;
    colnames.emplace_back("SHE_"+shearType+"_UPDATED_RA", "ra");
    colnames.emplace_back("SHE_"+shearType+"_UPDATED_DEC", "dec");
    colnames.emplace_back("SHE_"+shearType+"_G1", "g1");
    colnames.emplace_back("SHE_"+shearType+"_G2", "g2");
    colnames.emplace_back("SHE_"+shearType+"_WEIGHT", "w");
    colnames.emplace_back("PHZ_"+shearType+"_CORRECTION", "z_corr");
    colnames.emplace_back("PHZ_MEDIAN", "z");
    return colnames;
}

std::vector<std::pair<std::string, std::string>> getClusterColNamesAndProxy()
{
    std::vector<std::pair<std::string, std::string>> colnames;
    colnames.emplace_back("ID_DET_CLUSTER", "id");
    colnames.emplace_back("RIGHT_ASCENSION_CLUSTER", "ra");
    colnames.emplace_back("DECLINATION_CLUSTER", "dec");
    colnames.emplace_back("Z_CLUSTER", "z");
    colnames.emplace_back("RICHNESS_CLUSTER", "richness");
    colnames.emplace_back("RADIUS_CLUSTER", "r");
    return colnames;
}

int getIndexCol(std::vector<std::string>& colname, std::string string)
{
    int ind = -1;
    auto it = std::find(colname.begin(), colname.end(), string);
    if (it == colname.end())
    {
        logger.info() << "not in colnames";
    }
    else
    {
        ind = std::distance(colname.begin(), it);
        logger.info() << "column name " << string << " index is: " << ind;
    }
    return ind;
}

std::string getColName(std::vector<std::string>& colname,
        std::string substring1, std::string substring2)
{
    std::string name;
    for (size_t i = 0; i < colname.size(); i++)
    {
        if (colname[i].find(substring1) != std::string::npos
                || colname[i].find(substring2) != std::string::npos)
        {
            name = colname[i];
        }
    }
    return name;
}

std::string getColName(std::vector<std::string>& colname, std::string substring)
{
    std::string name;
    for (size_t i = 0; i < colname.size(); i++)
    {
        if (colname[i].find(substring) != std::string::npos)
        {
            name = colname[i];
        }
    }
    return name;
}

std::pair<Healpix_Map<double>, Healpix_Map<double> > readHealpixMap(
        const std::string &Map)
{
    arr<double> myarrE;
    arr<double> myarrB;

    MefFile fitsFile(Map, FileMode::Read);
    const auto &ext = fitsFile.access<BintableHdu>(1);
    const auto records = ext.header().parseAll();
    const auto Ttype1 = records.as<std::string>("TTYPE1");
    const auto Ttype2 = records.as<std::string>("TTYPE2");
    std::string ordering = records.as<std::string>("ORDERING");
    if (((Ttype1.value).compare("GAMMA1") == 0)
            || ((Ttype1.value).compare("KAPPA_E") == 0)
            || ((Ttype1.value).compare("PIXEL") == 0))
    {
        auto mE = ext.readColumn<double>(ColumnKey(Ttype1));
        myarrE.copyFrom(mE.vector());
    }
    if (((Ttype2.value).compare("GAMMA2") == 0)
            || ((Ttype2.value).compare("KAPPA_B") == 0)
            || ((Ttype2.value).compare("WEIGHT") == 0))
    {
        const auto mB = ext.readColumn<double>(ColumnKey(Ttype2));
        myarrB.copyFrom(mB.vector());
    }

    Healpix_Map<double> mapE, mapB;
    mapE.Set(myarrE, ordering == "RING" ? RING : NEST);
    mapB.Set(myarrB, ordering == "RING" ? RING : NEST);

    return std::pair<Healpix_Map<double>, Healpix_Map<double> >(mapE, mapB);
}

std::pair<double, double> getIndex2RaDec(Healpix_Map<double>& map, int index)
{
    pointing ptg;
    ptg = map.pix2ang(index);
    //logger.info() << "dec: " << -( ptg.theta - M_PI*double(0.5) )/deg2rad;
    //logger.info() << "ra: " << (ptg.phi/deg2rad + rotPhi);
    double dec = (-(ptg.theta - M_PI * double(0.5)) / deg2rad);
    double ra = (ptg.phi / deg2rad + rotPhi);
    return std::pair<double, double>(ra, dec);
}

void fillJsonOutput(const fs::path& workdir, const fs::path& outFilename,
        std::vector<fs::path>& inFilenames)
{
    std::ofstream outfile((workdir / outFilename).string(), std::ios_base::app);
    size_t nFiles = inFilenames.size();
    outfile << "[";
    for (size_t i = 0; i < nFiles; i++)
    {
        outfile << inFilenames[i].filename();
        if (i < nFiles - 1)
        {
            outfile << ", ";
        }
    }
    outfile << "]";
    outfile.close();
}

void makeGaussianKernel(double* kernel, double sigmaX, double sigmaY, int sizeX,
        int sizeY)
{
    double sigmaX2 = sigmaX * sigmaX;
    double sigmaY2 = sigmaY * sigmaY;

    // Fill in the gaussian kernel
    double sum(0.);
    for (int j = 0; j != sizeY; ++j)
    {
        float y = j - (sizeY - 1.) / 2.;
        for (int i = 0; i != sizeX; ++i)
        {
            float x = i - (sizeX - 1.) / 2.;
            kernel[j * sizeY + i] = exp(
                    -(x * x / (2 * sigmaX2) + y * y / (2 * sigmaY2)));
            sum += kernel[j * sizeY + i];
        }
    }
    // Normalize the kernel
    for (int k = 0; k != sizeX * sizeY; ++k)
    {
        kernel[k] /= sum;
    }
}

void replaceSubStringInPath(fs::path& path, std::string search,
        std::string substitute)
{
    std::string str = path.string();
    boost::replace_all(str, search, substitute);
    path = fs::path(str);
}

VecColumn<double> generateColumn(const std::string &name,
                                 std::vector<double>& data)
{
    ColumnInfo<double> info(name, "");
    return makeColumn(info, std::move(data));
}

void fillVecColumnLinspace(VecColumn<double>& col,
                           size_t N, double start, double stop,
                           bool endpoint)
{
    double step;
    if(endpoint)
    {
        step = (stop - start) / (N - 1);
    }
    else
    {
        step = (stop - start) / N ;
    }
    for (size_t i = 0; i < N; i++)
    {
        col(i) = start + step * i;
    }
}

std::string getXmlProductType(const fs::path& filepath)
{
    boost::property_tree::ptree tree;
    read_xml(filepath.native(), tree);
    std::string product_type = tree.begin()->first;
    product_type = splitString(product_type, ":").back();
    return product_type;
}

} // LE3_2D_MASS_WL_UTILITIES
