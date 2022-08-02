/**
 * @file src/lib/GetMap.cpp
 * @date 05/13/19
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

#include "LE3_2D_MASS_WL_CARTESIAN/GenericMap.h"

static Elements::Logging logger = Elements::Logging::getLogger("GenericMap");

using namespace LE3_2D_MASS_WL_UTILITIES;

namespace LE3_2D_MASS_WL_CARTESIAN
{
GenericMap::GenericMap() : m_sizeXaxis(0), m_sizeYaxis(0), m_sizeZaxis(0),
                           m_nGalaxies(0), m_PixelSize(0.001)
{
}

GenericMap::GenericMap(int Xaxis, int Yaxis, int Zaxis, int nbGalaxies,
                       double pixelSize) : m_sizeXaxis(Xaxis),
                       m_sizeYaxis(Yaxis), m_sizeZaxis(Zaxis),
                       m_nGalaxies(nbGalaxies), m_PixelSize(pixelSize)
{
    initMatrices();
}

GenericMap::GenericMap(const fs::path& filename) : m_sizeXaxis(0),
        m_sizeYaxis(0), m_sizeZaxis(0), m_nGalaxies(0),
        m_PixelSize(0.001)
{
    fillFromFitsFile(filename.native());
}

GenericMap::GenericMap(const std::string& filename) : m_sizeXaxis(0),
                       m_sizeYaxis(0), m_sizeZaxis(0), m_nGalaxies(0),
                       m_PixelSize(0.001)
{
    fillFromFitsFile(filename); // will init matrices size and read its header
}

GenericMap::GenericMap(GenericMap const& copyMap, bool copyValues) :
                       m_sizeXaxis(copyMap.m_sizeXaxis),
                       m_sizeYaxis(copyMap.m_sizeYaxis),
                       m_sizeZaxis(copyMap.m_sizeZaxis),
                       m_nGalaxies(copyMap.m_nGalaxies),
                       m_PixelSize(copyMap.m_PixelSize)
{
    initMatrices();

    // Assign values from the input array to the map
    if (copyValues)
    {
        fullCopy(copyMap);
    }
}

void GenericMap::fullCopy(GenericMap const& copyMap)
{
    for (int k = 0; k != m_sizeZaxis; ++k)
    {
        singleAxisCopy(copyMap, k);
    }
}

void GenericMap::singleAxisCopy(GenericMap const& copyMap, int k)
{
    for (int i = 0; i != m_sizeXaxis; ++i)
    {
        for (int j = 0; j != m_sizeYaxis; ++j)
        {
            m_mapValues[k](i, j) = copyMap.m_mapValues[k](i, j);
        }
    }
}

void GenericMap::applyThreshold(double threshold, int k)
{
    m_mapValues[k].applyThreshold(threshold);
}

double GenericMap::getSigma(int k) const
{
    return m_mapValues[k].getSigma();
}

double GenericMap::getMin(int k) const
{
    return m_mapValues[k].getMin();
}

double GenericMap::getMax(int k) const
{
    return m_mapValues[k].getMax();
}

double GenericMap::getFlux(int k) const
{
    return m_mapValues[k].getFlux();
}

void GenericMap::initMatrices()
{
    m_mapValues.resize(m_sizeZaxis);
    for (int i = 0; i < m_sizeZaxis; i++)
    {
        m_mapValues[i].resize(m_sizeXaxis, m_sizeYaxis);
        m_mapValues[i].clear();
    }
}

void GenericMap::updateSizes(int sizeXaxis, int sizeYaxis, int sizeZaxis)
{
    m_sizeXaxis = sizeXaxis;
    m_sizeYaxis = sizeYaxis;
    m_sizeZaxis = sizeZaxis;
    initMatrices();
}

void GenericMap::fillFromFitsFile(const std::string& filename, int hduIndex)
{
    if (!fs::exists(filename))
    {
        logger.fatal() << "File does not exist, map is not filled: "
                       << filename;
        return;
    }

    MefFile fitsfile(filename, FileMode::Edit);

    // List of HDU
    std::vector<std::string> Hdu_names = fitsfile.readHduNames();

    // Open hdu and fill
    const auto &ext = fitsfile.accessFirst<ImageHdu>(Hdu_names[hduIndex]);
    const auto ramin = ext.header().parseOr<double>({"RAMIN", 0}) * M_PI / 180.;
    const auto ramax = ext.header().parseOr<double>({"RAMAX", 0}) * M_PI / 180.;
    const auto decmin = ext.header().parseOr<double>({"DECMIN", 0}) * M_PI / 180.;
    const auto decmax = ext.header().parseOr<double>({"DECMAX", 0}) * M_PI / 180.;
    const auto zmin = ext.header().parseOr<double>({"ZMIN", 0});
    const auto zmax = ext.header().parseOr<double>({"ZMAX", 0});
    // m_CB = PatchDef(ramin, ramax, decmin, decmax, zmin, zmax);
    m_nGalaxies = ext.header().parseOr<int>({"NGAL", 0});
    m_PixelSize = ext.header().parseOr<double>({"PIXSIZE", 0}) * M_PI / 180.;

    logger.info() << "Pixel Size found in input Map (radians): " << m_PixelSize;
    logger.info() << "Pixel Size found in input Map (arcmin): "
                  << m_PixelSize * 180. / M_PI * 60;

    //read Image
    const auto image = ext.readRaster<double, 3>();
    const auto sizes = image.shape().vector();
    m_sizeXaxis = sizes[0];
    m_sizeYaxis = sizes[1];
    m_sizeZaxis = sizes[2];

    initMatrices();

    // Assign values from the input array to the map
    for (int z = 0; z != m_sizeZaxis; ++z)
    {
        for (int y = 0; y != m_sizeYaxis; ++y)
        {
            for (int x = 0; x != m_sizeXaxis; ++x)
            {
                m_mapValues[z](x, y) = image[{ x, y, z }];
            }
        }
    }
}

std::vector<double> GenericMap::getMeanValues() const
{
    // Create the vector of mean values
    std::vector<double> meanValues;

    // Loop over all the values of the map
    for (int k = 0; k != m_sizeZaxis; ++k)
    {
        double mean(0.);
        for (int i = 0; i != m_sizeXaxis; ++i)
        {
            for (int j = 0; j != m_sizeYaxis; ++j)
            {
                mean += m_mapValues[k](i, j);
            }
        }
        // Scale the sum to the mean
        mean /= m_sizeXaxis * m_sizeYaxis;

        // Add that value to the output vector
        meanValues.push_back(mean);
    }
    return meanValues;
}

void GenericMap::removeOffset(std::vector<double> offset)
{
    // For each value in the map, add the offset value
    for (int k = 0; k != m_sizeZaxis; ++k)
    {
        for (int i = 0; i != m_sizeXaxis; ++i)
        {
            for (int j = 0; j != m_sizeYaxis; ++j)
            {
                m_mapValues[k](i, j) -= offset[k];
            }
        }
    }
}

void GenericMap::wrapIndex(int& binx, int& biny, int& binz) const
{
    // Make sure the bin values are not wrong
    if (binx >= m_sizeXaxis)
    {
        binx = m_sizeXaxis - 1;
    }
    if (biny >= m_sizeYaxis)
    {
        biny = m_sizeYaxis - 1;
    }
    if (binz >= m_sizeZaxis)
    {
        binz = m_sizeZaxis - 1;
    }
}

unsigned int GenericMap::getNumberOfGalaxies() const
{
    return m_nGalaxies;
}

void GenericMap::setNumberOfGalaxies(unsigned int N)
{
    m_nGalaxies = N;
}

double GenericMap::getPixelSize() const
{
    return m_PixelSize;
}

double GenericMap::getBinValue(int binx, int biny, int binz) const
{
    wrapIndex(binx, biny, binz);
    return m_mapValues[binz](binx, biny);
}

void GenericMap::setBinValue(int binx, int biny, int binz, double value)
{
    wrapIndex(binx, biny, binz);
    m_mapValues[binz](binx, biny) = value;
}

Matrix& GenericMap::operator[](int k)
{
    return m_mapValues[k];
}

const Matrix& GenericMap::operator[](int k) const
{
    return m_mapValues[k];
}

double* GenericMap::getImageAddress(int k)
{
    return m_mapValues[k].data().begin();
}

void GenericMap::clear(int k)
{
    m_mapValues[k].clear();
}

void GenericMap::add_borders()
{
    // Backup the original array
    std::vector<Matrix> tmpArray(m_mapValues);

    // Update size
    m_sizeXaxis *= 2;
    m_sizeYaxis *= 2;
    for (int i = 0; i < m_sizeZaxis; i++)
    {
        m_mapValues[i].resize(m_sizeXaxis, m_sizeYaxis);
    }

    // Assign values from the temporary array to the bordered one
    for (int k = 0; k != m_sizeZaxis; ++k)
    {
        for (int i = 0; i != m_sizeXaxis; ++i)
        {
            for (int j = 0; j != m_sizeYaxis; ++j)
            {
                if (i < m_sizeXaxis / 4 || i >= 3 * m_sizeXaxis / 4
                 || j < m_sizeYaxis / 4 || j >= 3 * m_sizeYaxis / 4)
                {
                    m_mapValues[k](i, j) = 0.;
                }
                else
                {
                    m_mapValues[k](i, j) = tmpArray[k](i - m_sizeXaxis / 4,
                                                       j - m_sizeYaxis / 4);
                }
            }
        }
    }
}

void GenericMap::remove_borders()
{
    // Backup the bordered array
    std::vector<Matrix> tmpArray(m_mapValues);

    // Update size
    m_sizeXaxis /= 2;
    m_sizeYaxis /= 2;
    for (int i = 0; i < m_sizeZaxis; i++)
    {
        m_mapValues[i].resize(m_sizeXaxis, m_sizeYaxis);
    }

    // Assign values from the old array to the new without borders
    for (int k = 0; k != m_sizeZaxis; ++k)
    {
        for (int i = 0; i != m_sizeXaxis; ++i)
        {
            for (int j = 0; j != m_sizeYaxis; ++j)
            {
                m_mapValues[k](i, j) = tmpArray[k](i + m_sizeXaxis / 2,
                                                   j + m_sizeYaxis / 2);
            }
        }
    }
}

unsigned int GenericMap::getXdim() const
{
    return m_sizeXaxis;
}

unsigned int GenericMap::getYdim() const
{
    return m_sizeYaxis;
}

unsigned int GenericMap::getZdim() const
{
    return m_sizeZaxis;
}

VecRaster<double, 3> GenericMap::createRaster()
{
    // Create and fill a raster
    VecRaster<double, 3> raster(
    { m_sizeXaxis, m_sizeYaxis, m_sizeZaxis });
    for (int z = 0; z < m_sizeZaxis; ++z)
    {
        for (int y = 0; y < m_sizeYaxis; ++y)
        {
            for (int x = 0; x < m_sizeXaxis; ++x)
            {
                raster[{ x, y, z }] = m_mapValues[z](x, y);
            }
        }
    }
    return raster;
}

void GenericMap::setArray3D(double *array)
{
    for (int k = 0; k != m_sizeZaxis; ++k)
    {
        for (int i = 0; i != m_sizeXaxis; ++i)
        {
            for (int j = 0; j != m_sizeYaxis; ++j)
            {
                m_mapValues[k](i, j) = array[m_sizeXaxis * m_sizeYaxis * k
                                            + m_sizeYaxis * i + j];
            }
        }
    }
}

void GenericMap::getArray(double *array) const
{
    // Assign values to the elements
    for (int k = 0; k != m_sizeZaxis; ++k)
    {
        for (int i = 0; i != m_sizeXaxis; ++i)
        {
            for (int j = 0; j != m_sizeYaxis; ++j)
            {
                array[m_sizeXaxis * m_sizeYaxis * k + m_sizeYaxis * i + j] =
                        m_mapValues[k](i, j);
            }
        }
    }
}

void GenericMap::operator=(const GenericMap& map)
{
    m_sizeXaxis = map.m_sizeXaxis;
    m_sizeYaxis = map.m_sizeYaxis;
    m_sizeZaxis = map.m_sizeZaxis;
    m_nGalaxies = map.m_nGalaxies;
    m_PixelSize = map.m_PixelSize;
    m_mapValues = map.m_mapValues;
}

bool GenericMap::pixelate(int xBinning, int yBinning)
{
    if (xBinning == 0 && yBinning == 0)
    {
        return false;
    }

    // Put the binning from power of two to actual values
    xBinning = pow(2, xBinning);
    yBinning = pow(2, yBinning);

    // Check the asked binning is not higher or equal to the number of current bins
    if (xBinning >= m_sizeXaxis || yBinning >= m_sizeYaxis)
    {
        return false;
    }

    // Backup the array
    std::vector<Matrix> tmpArray(m_mapValues);

    // Resize and clear
    for (int i = 0; i < m_sizeZaxis; i++)
    {
        m_mapValues[i].resize(m_sizeXaxis / xBinning, m_sizeYaxis / yBinning);
        m_mapValues[i].clear();
    }

    // Reassign binned values to the map
    for (int k = 0; k != m_sizeZaxis; ++k)
    {
        for (int i = 0; i != m_sizeXaxis; ++i)
        {
            for (int j = 0; j != m_sizeYaxis; ++j)
            {
                m_mapValues[k](i / xBinning, j / yBinning) += tmpArray[k](i, j)
                        / xBinning / yBinning;
            }
        }
    }

    m_sizeXaxis /= xBinning;
    m_sizeYaxis /= yBinning;

    return true;
}

void GenericMap::writeImageHeader(const Hdu &hdu,
                                  CartesianParam &cartesianParam)
{
    auto& patch = cartesianParam.getPatches()[0];
    float cd11 = (patch.getRaMax()-patch.getRaMin())*180./M_PI/m_sizeXaxis;
    float cd22 = (patch.getDecMax()-patch.getDecMin())*180./M_PI/m_sizeYaxis;
    float cd12 = 0.0;
    float cd21 = 0.0;
    float crval1 = (floor(m_sizeXaxis/2)+0.5)*cd11+patch.getRaMin()*180./M_PI;
    float crval2 = (floor(m_sizeYaxis/2)+0.5)*cd22+patch.getRaMax()*180./M_PI;

    std::list<Record<boost::any>> records =
    {
    { "WCSAXES", 2, "", "Number of axes in World Coordinate System" },
    { "CRPIX1", "", "", "Pixel coordinate of reference point" },
    { "CRPIX2", "", "", "Pixel coordinate of reference point" },
    { "PC1_1", cd11, "", "Coordinate transformation matrix element" },
    { "PC1_2", cd12, "", "Coordinate transformation matrix element" },
    { "PC2_1", cd21, "", "Coordinate transformation matrix element" },
    { "PC2_2", cd22, "", "Coordinate transformation matrix element" },
    { "CDELT1", "", "deg", "Coordinate increment at reference point" },
    { "CDELT2", "", "deg", "Coordinate increment at reference point" },
    { "CUNIT1", "deg", "", "Unit of the first coordinate value" },
    { "CUNIT2", "deg", "", "Unit of the second coordinate value" },
    { "CTYPE1", "RA---TAN", "", "Right ascension, gnomonic projection" },
    { "CTYPE2", "DEC--TAN", "", "Declination, gnomonic projection" },
    { "CRVAL1", crval1, "deg", "Coordinate value at reference point" },
    { "CRVAL2", crval2, "deg", "Coordinate value at reference point" },
    { "LONPOLE", "", "deg", "Native longitude of celestial pole" },
    { "LATPOLE", "", "deg", "Native latitude of celestial pole" },
    { "RADESYS", "", "", "Equatorial coordinate system" },
    { "EQUINOX", "", "", "Equinox of celestial coordinate system (e.g. 2000)" },
    { "DATE-OBS", "", "", "Start of observation" }, //format yyyy-mm-ddThh:mm:ss.sss
    { "DATE-END", "", "", "End of observation" },
    { "SHE_SOFT", "LENSMC", "", "Origin of Shear Catalog" },
    { "SHE_SEL", "FITCLASS=0", "", "Origin of Shear Catalog" },
    { "SOFTNAME", "LE3_2D_MASS_WL_KS", "", "Software used to create the product" },
    { "SOFTVERS", SWVersion, "", "Software version" },
    { "PROJ", "TAN", "", "projection used" },
    { "PWIDTH", patch.getPatchWidth()*rad2deg, "deg", "PatchWidth" },
    { "PIXSIZE", patch.getPixelSize()*rad2deg, "deg", "Pixelsize" },
    { "NITREDSH", cartesianParam.getRsNItReducedShear(), "",
            "Number of iterations for reduced shear" },
    { "STDREDSH", cartesianParam.getRsGaussStd(), "",
            "gaussian smoothing sigma for reduced shear" },
    { "NITINP", cartesianParam.getNInpaint(), "",
            "Number of iterations for inpainting" },
    { "NSCINP", cartesianParam.getNInpScale(),
            "Number of scales for inpainting" },
    { "VARPERSC", cartesianParam.isEqualVarPerScale(),
            "True if equal variance per scale forced in" },
    { "FBMODE", cartesianParam.isForceBMode(), "",
            "True if B-mode forced to zero in the gaps" },
    { "ADDBORD", cartesianParam.isAddBorder(), "", "True if Borders added" },
    { "DENTYPE", "GAUSSIAN", "", "denoising type" },
    { "GAUSSSTD", cartesianParam.getGaussStd(), "",
            "Standard deviation of gaussian smoothing" },
    { "FDRVAL", cartesianParam.getThresholdFdr(), "",
            "false discovery rate threshold" },
    { "NZBINS", cartesianParam.getNbins(), "", "number of redshift bins" },
    { "BAL_BINS", cartesianParam.isBalancedBins(), "",
            "True if balanced bins are applied" },
    { "NRESAMPL", cartesianParam.getNResamples(), "",
            "number of resampling of input dictionary" },
    // TODO: fix this!
    //{ "ZMIN", patch.getZMin(), "", "Minimum Redshift" },
    //{ "ZMAX", patch.getZMax(), "", "Maximum Redshift" }
    };

    logger.info() << "Record created";

    hdu.header().writeSeq(records);

}

void GenericMap::writeMap(const std::string& filename,
                          CartesianParam &cartesianParam)
{
    MefFile f(filename, FileMode::Overwrite);
    const auto raster = createRaster();
    logger.info() << "writing records to primary HDU";
    const auto &primary = f.primary();
    writeImageHeader(primary, cartesianParam);
    logger.info() << "Assigning new Image HDU: " << cartesianParam.getExtName();
    const auto &ext = f.assignImageExt(cartesianParam.getExtName(), raster);

    // Image HDU
    //Adding extra keys to image header
    auto& patch = cartesianParam.getPatches()[0];
    std::list<Record<boost::any>> records =
    {
    //TODO:fix this! (need param)
    //{"ZMIN", patch.getZMin()},
    //{"ZMAX", patch.getZMax()},
    {"RAMIN", patch.getRaMin() * rad2deg},
    {"RAMAX", patch.getRaMax() * rad2deg},
    {"DECMIN", patch.getDecMin() * rad2deg},
    {"DECMAX", patch.getDecMax() * rad2deg},
    {"NGAL", int(getNumberOfGalaxies())},
    {"PIXSIZE", patch.getPixelSize()}
    };
    ext.header().writeSeq(records);
}

void GenericMap::writeMap(fs::path filename, CartesianParam &cartesianParam)
{
    writeMap(filename.native(), cartesianParam);
}

}  // namespace LE3_2D_MASS_WL_CARTESIAN
