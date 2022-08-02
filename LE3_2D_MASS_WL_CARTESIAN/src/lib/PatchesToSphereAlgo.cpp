/**
 * @file src/lib/PatchesToSphereAlgo.cpp
 * @date 08/17/21
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

#include "LE3_2D_MASS_WL_CARTESIAN/PatchesToSphereAlgo.h"

static Elements::Logging logger = Elements::Logging::getLogger(
        "PatchesToSphereAlgo");

using namespace LE3_2D_MASS_WL_UTILITIES;

namespace LE3_2D_MASS_WL_CARTESIAN
{

PatchesToSphereAlgo::PatchesToSphereAlgo(CartesianParam &params) :
        m_cartesianParam(params)
{
}

bool PatchesToSphereAlgo::getNoisyConvergenceMaps(fs::path& workdir,
                                                  fs::path& InShearMap,
                                                  fs::path& outputMaps)
{
    // Object to perform Cartesian KS algorithm
    CartesianAlgoKS CartesainAlgo(m_cartesianParam);
    float SigmaGauss = 0.;
    std::vector<fs::path> filenames = readFilenamesInJson(workdir / InShearMap);
    fs::path datadir = workdir / "data";
    if (m_cartesianParam.getGaussStd() != 0)
    {
        SigmaGauss = m_cartesianParam.getGaussStd();
    }
    m_cartesianParam.setGaussStd(0.);
    size_t pos;
    std::ofstream outfile;
    outfile.open((workdir / outputMaps).string(), std::ios_base::app);
    outfile << "[";
    for (size_t i = 0; i < filenames.size(); i++)
    {
        logger.info() << filenames[i];
        pos = (filenames[i].string()).find("ShearMap");
        fs::path outConvergenceMap(
                "EUC_LE3_WL_NoisyConvergenceMapKS_"
                        + (filenames[i].string()).substr(pos));

        ShearMap shearMap(datadir / filenames[i]);
        ConvergenceMap convergenceMap(shearMap, false);
        CartesainAlgo.performMassMapping(shearMap, convergenceMap, workdir);
        convergenceMap.writeMap(datadir / outConvergenceMap, m_cartesianParam);

        outfile << outConvergenceMap.filename();
        if (i < filenames.size() - 1)
        {
            outfile << ",";
        }
    } //end of for loop
    outfile << "]";
    outfile.close();
    m_cartesianParam.setGaussStd(SigmaGauss);
    return true;
}

bool PatchesToSphereAlgo::getDeNoisyConvergenceMaps(fs::path& workdir,
                                                    fs::path& InShearMap,
                                                    fs::path& outputMaps)
{
    // Object to perform Cartesian KS algorithm
    CartesianAlgoKS CartesainAlgo(m_cartesianParam);
    fs::path outConvergenceMap
    { };
    std::vector<fs::path> filenames = readFilenamesInJson(workdir / InShearMap);
    fs::path datadir = workdir / "data";
    size_t pos;
    std::ofstream outfile;
    outfile.open((workdir / outputMaps).string(), std::ios_base::app);
    outfile << "[";

    for (size_t i = 0; i < filenames.size(); i++)
    {
        pos = (filenames[i].string()).find("ShearMap");
        outConvergenceMap = fs::path(
                "EUC_LE3_WL_DenoisedConvergenceMapKS_"
                        + (filenames[i].string()).substr(pos));

        ShearMap shearMap(datadir / filenames[i]);
        ConvergenceMap convergenceMap(shearMap, false);
        CartesainAlgo.performMassMapping(shearMap, convergenceMap, workdir);

        outfile << outConvergenceMap.filename();
        if (i < filenames.size() - 1)
        {
            outfile << ",";
        }
    } //end of for loop

    outfile << "]";
    outfile.close();
    return true;
}

void PatchesToSphereAlgo::getPatchData(fs::path &convMaps, fs::path &workdir,
                              const std::string& outCenterPatchPositionFileFits,
                              std::vector<std::vector<double>>& data)
{
    std::vector<double> ra_ctr;
    std::vector<double> dec_ctr;
    std::vector<double> ra;
    std::vector<double> dec;
    std::vector<double> PixelDataE;
    std::vector<double> PixelDataB;
    Projection projection;
    fs::path datadir = workdir / "data";
    std::vector<fs::path> filenames = readFilenamesInJson(workdir / convMaps);
    for (size_t it = 0; it < filenames.size(); it++)
    {
        logger.info() << "Reading Patch: " << datadir / filenames[it];
        ConvergenceMap convMap(datadir / filenames[it]);
        //TODO: fill patch from info in map fits
        PatchDef patch(0,0,0,0);
        double raMin = patch.getRaMin();
        double raMax = patch.getRaMax();
        double decMin = patch.getDecMin();
        double decMax = patch.getDecMax();
        unsigned int x = convMap.getXdim();
        unsigned int y = convMap.getYdim();

        double ra0 = 0.5 * (raMin + raMax); // in degree (Center of patch)
        double dec0 = 0.5 * (decMin + decMax); // in degree (Center of patch)

        ra_ctr.push_back(ra0);
        dec_ctr.push_back(dec0);
        logger.info() << "getCenterPatchPositionsFile";
        getCenterPatchPositionsFile(outCenterPatchPositionFileFits, ra_ctr,
                dec_ctr, patch);

        double raRange = (raMax - raMin) * M_PI / 180.; //in rad
        double decRange = (decMax - decMin) * M_PI / 180.; //in rad

        // Loop over all pixels
        for (unsigned int i = 0; i < x; i++)
        {
            for (unsigned int j = 0; j < y; j++)
            {
                // Perform transform from pixel location to ra and dec
                double tmpx = (i + 0.5) * raRange / x - 0.5 * raRange;
                double tmpy = (j + 0.5) * decRange / y - 0.5 * decRange;
                std::pair<double, double> radec =
                        projection.getInverseGnomonicProjection(tmpx, tmpy, ra0,
                                dec0);
                ra.push_back(radec.first);
                dec.push_back(radec.second);
                PixelDataE.push_back(convMap.getBinValue(i, j, 0));
                PixelDataB.push_back(convMap.getBinValue(i, j, 1));
            }
        }
    }
    data.push_back(ra);
    data.push_back(dec);
    data.push_back(PixelDataE);
    data.push_back(PixelDataB);
}

void PatchesToSphereAlgo::getHealpixFormatMap(
        std::vector<std::vector<double>>& data, Healpix_Map<double>& mapE,
        Healpix_Map<double>& mapB, Healpix_Map<double>& mapGalCount)
{
    int order = hb.nside2order(m_cartesianParam.getNside());
    mapE.Set(order, RING);
    mapB.Set(order, RING);
    mapGalCount.Set(order, RING);

    mapE.fill(0.);
    mapB.fill(0.);
    mapGalCount.fill(0.);

    unsigned int npix = mapE.Npix();
    logger.info() << "npix: " << npix;
    int ngal = data[0].size();
    logger.info() << "ngal: " << ngal;

    for (int i = 0; i < ngal; i++)
    {
        double theta = -M_PI / 180.0 * data[1][i] + M_PI * double(0.5);
        double phi = M_PI / 180.0 * (data[0][i]);
        pointing ptg = pointing(theta, phi);
        ptg.normalize();
        auto id = mapE.ang2pix(ptg);
        mapE[id] += data[2][i];
        mapB[id] += data[3][i];
        mapGalCount[id] += 1.;
    }
    for (unsigned int id_pix = 0; id_pix < npix; id_pix++)
    {
        if (mapGalCount[id_pix] != 0)
        {
            mapE[id_pix] /= mapGalCount[id_pix];
            mapB[id_pix] /= mapGalCount[id_pix];
        }
    }
}

bool PatchesToSphereAlgo::getCenterPatchPositionsFile(
        const std::string& outFileFits, std::vector<double>& ra_ctr,
        std::vector<double>& dec_ctr, PatchDef& m_CB)
{
    MefFile f(outFileFits, FileMode::Overwrite);

    logger.info() << "writing fits file";
    const auto col1 = generateColumn("RA_CTR", ra_ctr);
    const auto col2 = generateColumn("DEC_CTR", dec_ctr);

    const auto &ext = f.assignBintableExt("PATCHES_POSITIONS", col1, col2);
    //Adding extra keys to image header
    std::list<Record<double>> records =
    {
    //TODO: fix this!
    //{ "ZMIN", m_CB.getZMin() },
    //{ "ZMAX", m_CB.getZMax() }
    };
    ext.header().writeSeq(records);

    return true;
}

void PatchesToSphereAlgo::writePrimaryHeader(const Hdu &hdu)
{
    std::list<Record<boost::any>> records =
    {
    { "RADESYS", "", "", "Equatorial coordinate system" },
    { "EQUINOX", "", "", "Equinox of celestial coordinate system (e.g. 2000)" },
    { "DATE-OBS", "", "", "Start of observation" }, //format yyyy-mm-ddThh:mm:ss.sss
            { "DATE-END", "", "", "End of observation" },
            { "SHE_SOFT", "LENSMC", "", "Origin of Shear Catalog" },
            { "SHE_SEL", "FITCLASS=0", "", "Origin of Shear Catalog" },
            { "SOFTNAME", "LE3_2D_MASS_WL_KS", "",
                    "Software used to create the product" },
            { "SOFTVERS", SWVersion, "", "Software version" },
            { "NITREDSH", m_cartesianParam.getRsNItReducedShear(), "",
                    "Number of iterations for reduced shear" },
            { "NITINP", m_cartesianParam.getNInpaint(), "",
                    "Number of iterations for inpainting" },
            { "NSCINP", m_cartesianParam.getNInpScale(),
                    "Number of scales for inpainting" },
            { "VARPERSC", m_cartesianParam.isEqualVarPerScale(),
                    "True if equal variance per scale forced in" },
            { "FBMODE", m_cartesianParam.isForceBMode(), "",
                    "True if B-mode forced to zero in the gaps" },
            { "DENTYPE", "GAUSSIAN", "", "denoising type" },
            { "GAUSSSTD", m_cartesianParam.getGaussStd(), "",
                    "Standard deviation of gaussian smoothing" },
            { "FDRVAL", m_cartesianParam.getThresholdFdr(), "",
                    "false discovery rate threshold" },
            { "NZBINS", m_cartesianParam.getNbins(), "",
                    "number of redshift bins" },
            { "BAL_BINS", m_cartesianParam.isBalancedBins(), "",
                    "True if balanced bins are applied" },
            { "NRESAMPL", m_cartesianParam.getNResamples(), "",
                    "number of resampling of input dictionary" },
            { "MCREAL", "", "", "" },
            { "MCSEED", "", "", "" }, };
    hdu.header().writeSeq(records);
}

void PatchesToSphereAlgo::writeHdu(const Hdu &hdu, int Nside)
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

bool PatchesToSphereAlgo::write_Map(const std::string& filename,
        Healpix_Map<double>& map, const std::string &colname)
{
    arr<double> myarr = map.Map();
    std::vector<double> mapValues;
    myarr.copyTo(mapValues);
    if (colname.compare("GALCOUNT") == 0)
    {
        std::string GalExt = "GALCOUNT_PATCHES";
        m_cartesianParam.setExtName(GalExt);
    }
    mapValues.shrink_to_fit();
    VecColumn<double> col = generateColumn(colname, mapValues);

    if (access(filename.c_str(), F_OK) != -1)
    {
        MefFile P2SFile(filename, FileMode::Edit);
        long nbHdu = P2SFile.hduCount();
        const auto &ext = P2SFile.access<BintableHdu>(nbHdu - 1);
        const auto &btCol = ext.columns();
        btCol.init(col.info());
        btCol.write(col);
    }
    else
    {
        MefFile P2SFile(filename, FileMode::Overwrite);
        const auto &primary = P2SFile.primary();
        writePrimaryHeader(primary);
        const auto &ext = P2SFile.assignBintableExt(
                m_cartesianParam.getExtName(), col);
        writeHdu(ext, map.Nside());
    }

    return true;
}

} // namespace LE3_2D_MASS_WL_CARTESIAN

