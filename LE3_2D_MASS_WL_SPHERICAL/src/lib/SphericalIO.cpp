/**
 * @file src/lib/SphericalIO.cpp
 * @date 12/24/20
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

#include "LE3_2D_MASS_WL_SPHERICAL/SphericalIO.h"

using namespace Euclid::WeakLensing::TwoDMass;
using namespace LE3_2D_MASS_WL_UTILITIES;

namespace fs = boost::filesystem;

static Elements::Logging logger = Elements::Logging::getLogger("SphericalIO");

namespace Euclid
{
namespace WeakLensing
{
namespace TwoDMass
{
namespace Spherical
{

SphericalIO::SphericalIO(LE3_2D_MASS_WL_SPHERICAL::SphericalParam &SphParam) :
        m_SphParam(SphParam)
{
}

void SphericalIO::writePrimaryHeader(const Hdu &hdu)
{
    std::list<Record<boost::any>> records =
    {
    { "RADESYS", "", "", "Equatorial coordinate system" },
    { "EQUINOX", "", "", "Equinox of celestial coordinate system (e.g. 2000)" },
    { "DATE-OBS", "", "", "Start of observation" },
    { "DATE-END", "", "", "End of observation" },
    { "SHE_SOFT", "LENSMC", "", "Origin of Shear Catalog" },
    { "SHE_SEL", "FITCLASS=0", "", "Origin of Shear Catalog" },
    { "SOFTNAME", "LE3_2D_MASS_WL_KS", "",
             "Software used to create the product" },
    { "SOFTVERS", SWVersion, "", "Software version" },
    { "NITREDSH", m_SphParam.getRsNItReducedShear(), "",
             "Number of iterations for reduced shear" },
    { "STDREDSH", m_SphParam.getRsGaussStd(), "",
             "gaussian smoothing sigma for reduced shear" },
    { "NITINP", m_SphParam.getNInpaint(), "",
             "Number of iterations for inpainting" },
    { "NSCINP", m_SphParam.getNInpScale(),
            "Number of scales for inpainting" },
    { "VARPERSC", m_SphParam.isEqualVarPerScale(),
            "True if equal variance per scale forced in" },
    { "FBMODE", m_SphParam.isForceBMode(), "",
             "True if B-mode forced to zero in the gaps" },
    { "DENTYPE", "GAUSSIAN", "", "denoising type" },
    { "GAUSSSTD", m_SphParam.getGaussStd(), "",
             "Standard deviation of gaussian smoothing" },
    { "FDRVAL", m_SphParam.getThresholdFdr(), "",
             "false discovery rate threshold" },
    { "NZBINS", m_SphParam.getNbins(), "", "number of redshift bins" },
    { "BAL_BINS", m_SphParam.isBalancedBins(), "",
             "True if balanced bins are applied" },
    { "NRESAMPL", m_SphParam.getNResamples(), "",
             "number of resampling of input dictionary" },
    { "MCREAL", "", "", "" },
    { "MCSEED", "", "", "" }, };
    hdu.header().writeSeq(records);
}

void SphericalIO::writeHdu(const Hdu &hdu, int Nside)
{
    std::list<Record<boost::any>> records =
    {
    { "PIXTYPE", "HEALPIX", "", "HEALPIX Pixelisation" },
    { "ORDERING", "RING", "", "Pixel ordering scheme, either RING or NESTED" },
    { "NSIDE", Nside, "", "Resolution parameter of HEALPIX" },
    { "FIRSTPIX", "0", "", "First pixel # (0 based)" },
    { "LASTPIX", "", "", "Last pixel # (0 based)" },
    { "INDXSCHM", "", "", "Indexing: IMPLICIT or EXPLICIT" },
    { "OBJECT", "0", "", "Sky coverage, either FULLSKY or PARTIAL" },
    { "ZMIN", "", "", "" },
    { "ZMAX", "", "", "" }, };
    hdu.header().writeSeq(records);
}

bool SphericalIO::write_Map(const std::string& filename,
        Healpix_Map<double>& map, const std::string &colname)
{
    arr<double> myarr = map.Map();
    std::vector<double> mapValues;
    myarr.copyTo(mapValues);
    if (colname.compare("GALCOUNT") == 0)
    {
        std::string GalExt = "GALCOUNT_SPHERE";
        m_SphParam.setExtName(GalExt);
    }
    //mapValues.erase( std::remove(mapValues.begin(), mapValues.end(), 0), mapValues.end()); //removing Zero values
    mapValues.shrink_to_fit();

    VecColumn<double> col = generateColumn(colname, mapValues);

    if (access(filename.c_str(), F_OK) != -1)
    {
        MefFile SphFile(filename, FileMode::Edit);
        long nbHdu = SphFile.hduCount();
        int nh = nbHdu - 1;
        const auto &ext = SphFile.access<BintableHdu>(nh);
        const auto &btCol = ext.columns();
        btCol.init(col.info());
        btCol.write(col);
    }
    else
    {
        MefFile SphFile(filename, FileMode::Overwrite);
        const auto &primary = SphFile.primary();
        writePrimaryHeader(primary);
        const auto &ext = SphFile.assignBintableExt(m_SphParam.getExtName(),
                col);
        writeHdu(ext, map.Nside());
    }

    return true;
}

bool SphericalIO::writeSphericalXMLfile(boost::filesystem::path& workdir,
        boost::filesystem::path& outConvergenceMapJson,
        boost::filesystem::path& GalCountJson,
        boost::filesystem::path& outXMLConvergenceMap)
{
    DmOutput dm;

    fs::path file = workdir / outXMLConvergenceMap;

    std::vector<fs::path> filenames = readFilenamesInJson(
            workdir / outConvergenceMapJson);

    std::string product_type = "DpdTwoDMassConvergenceSphere";
    auto product =
            initProduct<
                    dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere,
                    pro::le3::wl::twodmass::twoDMassCollectConvergenceSphere,
                    int>(product_type, m_SphParam.getNResamples());
    std::vector<fs::path> GalCountfilename = readFilenamesInJson(
            workdir / GalCountJson);
    fs::path galMapOut(GalCountfilename[0].native());
    dm.createSphereGalCountXml(product, galMapOut);
    for (size_t i = 0; i < filenames.size(); i++)
    {
        fs::path mapOut(filenames[i].native());
        if ((filenames[i].string()).find("NoisySphereConvergence")
                != std::string::npos)
        {
            dm.createNoisedSphereXml(product, mapOut);
        }
        if ((filenames[i].string()).find("DenoisedSphereConvergence")
                != std::string::npos)
        {
            dm.createDenoisedSphereXml(product, mapOut);
        }
    }
    writeProduct<
            dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere>(
            product, file);
    logger.info() << "DM output products created in: " << file;
    return true;
}

bool SphericalIO::writeSphericalMCXMLfile(boost::filesystem::path& workdir,
        boost::filesystem::path& outConvergenceMapJson,
        boost::filesystem::path& outXMLConvergenceMap)
{
    DmOutput dm;

    fs::path file = workdir / outXMLConvergenceMap;

    std::vector<fs::path> filenames = readFilenamesInJson(
            workdir / outConvergenceMapJson);

    std::string product_type = "DpdTwoDMassConvergenceSphere";
    auto product =
            initProduct<
                    dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere,
                    pro::le3::wl::twodmass::twoDMassCollectConvergenceSphere,
                    int>(product_type, m_SphParam.getNResamples());

    for (size_t i = 0; i < filenames.size(); i++)
    {
        fs::path mapOut(filenames[i].native());
        if ((filenames[i].string()).find("MCSphereConvergence")
                != std::string::npos)
        {
            dm.createSphereMCXml(product, mapOut);
        }
    }
    writeProduct<
            dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere>(
            product, file);
    logger.info() << "DM output products created in: " << file;
    return true;
}

} /* namespace Spherical */
} /* namespace TwoDMass */
} /* namespace WeakLensing */
} /* namespace Euclid */
