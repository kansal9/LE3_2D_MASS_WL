/**
 * @file LE3_2D_MASS_WL_UTILITIES/GenericParam.h
 * @date 02/04/22
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

#ifndef _LE3_2D_MASS_WL_UTILITIES_GENERICPARAM_H
#define _LE3_2D_MASS_WL_UTILITIES_GENERICPARAM_H

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_UTILITIES/PatchDef.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"

namespace LE3_2D_MASS_WL_UTILITIES {

/**
 * @class GenericParam
 * @brief
 *
 */

class  GenericParam {

public:
    /**
    * @brief Constructor
    */
    GenericParam();

    /**
     * @brief print all parameters
     */
    void print() const;

    /**
     * @brief get addBorder
     * @return addBorder
     */
    bool isAddBorder() const;

    /**
     * @brief set addBorder
     * @param addBorder
     */
    void setAddBorder(bool addBorder);

    /**
     * @brief get balancedBins
     * @return balancedBins
     */
    bool isBalancedBins() const;

    /**
     * @brief set balancedBins
     * @param balancedBins
     */
    void setBalancedBins(bool balancedBins);

    /**
     * @brief get decMax
     * @param decMin
     * @return decMax
     */
    double getDecMax(double decMin);

    /**
     * @brief set decMax
     * @param decMax
     */
    void setDecMax(double decMax);

    /**
     * @brief get decMin
     * @param mapCenterY
     * @return decMin
     */
    double getDecMin(double mapCenterY);

    /**
     * @brief set decMin
     * @param decMin
     */
    void setDecMin(double decMin);

    /**
     * @brief get denoisingAlgo
     * @return
     */
    const std::string& getDenoisingAlgo() const;

    /**
     * @brief get denoisingAlgo
     * @param denoisingAlgo
     */
    void setDenoisingAlgo(const std::string& denoisingAlgo);

    /**
     * @brief get equalVarPerScale
     * @return equalVarPerScale
     */
    bool isEqualVarPerScale() const;

    /**
     * @brief set equalVarPerScale
     * @param equalVarPerScale
     */
    void setEqualVarPerScale(bool equalVarPerScale);

    /**
     * @brief get extName
     * @return extName
     */
    const std::string& getExtName() const;

    /**
     * @brief set extName
     * @param extName
     */
    void setExtName(const std::string& extName);

    /**
     * @brief get forceBMode
     * @return forceBMode
     */
    bool isForceBMode() const;

    /**
     * @brief set forceBMode
     * @param forceBMode
     */
    void setForceBMode(bool forceBMode);

    /**
     * @brief get gaussStd
     * @return gaussStd
     */
    double getGaussStd() const;

    /**
     * @brief set gaussStd
     * @param gaussStd
     */
    void setGaussStd(double gaussStd);

    /**
     * @brief get massThreshold
     * @return massThreshold
     */
    double getMassThreshold() const;

    /**
     * @brief set massThreshold
     * @param massThreshold
     */
    void setMassThreshold(double massThreshold);

    /**
     * @brief get nbins
     * @return nbins
     */
    int getNbins() const;

    /**
     * @brief set nbins
     * @param nbins
     */
    void setNbins(int nbins);

    /**
     * @brief get nInpaint
     * @return nInpaint
     */
    int getNInpaint() const;

    /**
     * @brief set nInpaint
     * @param nInpaint
     */
    void setNInpaint(int nInpaint);

    /**
     * @brief get nInpScale
     * @return nInpScale
     */
    int getNInpScale() const;

    /**
     * @brief set nInpScale
     * @param nInpScale
     */
    void setNInpScale(int nInpScale);

    /**
     * @brief get nPatches
     * @return nPatches
     */
    unsigned int getNPatches() const;

    /**
     * @brief set nPatches
     * @param nPatches
     */
    void setNPatches(int nPatches);

    /**
     * @brief get nResamples
     * @return nResamples
     */
    int getNResamples() const;

    /**
     * @brief set nResamples
     * @param nResamples
     */
    void setNResamples(int nResamples);

    /**
     * @brief get nside
     * @return nside
     */
    int getNside() const;

    /**
     * @brief set nside
     * @param nside
     */
    void setNside(int nside);

    /**
     * @brief get overlap
     * @return overlap
     */
    int getOverlap() const;

    /**
     * @brief set overlap
     * @param overlap
     */
    void setOverlap(int overlap);

    /**
     * @brief get getParaFileType
     * @return getParaFileType
     */
    const std::string& getParaFileType() const;

    /**
     * @brief set getParaFileType
     * @param paraFileType
     */
    void setParaFileType(const std::string& paraFileType);

    /**
     * @brief get project
     * @return project
     */
    const std::string& getProject() const;

    /**
     * @brief set project
     * @param project
     */
    void setProject(const std::string& project);

    /**
     * @brief get rsGaussStd
     * @return rsGaussStd
     */
    double getRsGaussStd() const;

    /**
     * @brief set rsGaussStd
     * @param rsGaussStd
     */
    void setRsGaussStd(double rsGaussStd);

    /**
     * @brief get rsnItReducedShear
     * @return rsnItReducedShear
     */
    unsigned int getRsNItReducedShear() const;

    /**
     * @brief set rsnItReducedShear
     * @param rsnItReducedShear
     */
    void setRsNItReducedShear(int rsnItReducedShear);

    /**
     * @brief get rsThreshold
     * @return rsThreshold
     */
    double getRsThreshold() const;

    /**
     * @brief set rsThreshold
     * @param rsThreshold
     */
    void setRsThreshold(double rsThreshold);

    /**
     * @brief get thresholdFdr
     * @return thresholdFdr
     */
    double getThresholdFdr() const;

    /**
     * @brief set thresholdFdr
     * @param thresholdFdr
     */
    void setThresholdFdr(double thresholdFdr);

    /**
     * @brief get zMargin
     * @return zMargin
     */
    double getZMargin() const;

    /**
     * @brief set zMargin
     * @param zMargin
     */
    void setZMargin(double zMargin);

    /**
     * @brief get zMax
     * @return zMax
     */
    double getZMax(int i) const;

    /**
     * @brief get zMaxHalo
     * @return zMaxHalo
     */
    double getZMaxHalo() const;

    /**
     * @brief set zMaxHalo
     * @param zMaxHalo
     */
    void setZMaxHalo(double zMaxHalo);

    /**
     * @brief get zMin
     * @return zMin
     */
    double getZMin(int i) const;

    /**
     * @brief get Patches
     * @return Patches
     */
    const std::vector<PatchDef>& getPatches() const;

    /**
     * @brief get Patch
     * @return Patch
     */
    const PatchDef& getPatch(int i) const;

    /**
     * @brief readBaseParams
     * @param data_xml
     */
    template<typename T>
    void readBaseParams(T data_xml);

    /**
     * @brief readRedshiftBaseParams
     * @param data_xml
     */
    template<typename T>
    void readRedshiftBaseParams(T data_xml, CatalogData& cat);

protected:
    /**
     * number of iterations used to compute the reduced shear
     */
    unsigned int m_RSNItReducedShear;

    /**
     * threshold to perform hard-thresholding in the reduced shear
     */
    double m_RSThreshold;

    /**
     * standard deviation for the Gaussian filter (for reduced shear)
     */
    double m_RSGaussSTD;

    /**
     * name of denoising algo [GaussFilter]
     */
    std::string m_DenoisingAlgo;

    /**
     * standard deviation for the Gaussian filter
     */
    double m_GaussSTD;

    /**
     * False Discovery Rate (FDR) threshold if applied
     */
    double m_ThresholdFDR;

    /**
     * number of resampling that needs to be performed for Noise/SNR maps
     */
    unsigned int m_NResamples;

    /**
     * min redshifts for the patches
     */
    std::vector<double> m_ZMin;

    /**
     * max redshifts for the patches
     */
    std::vector<double> m_ZMax;

    /**
     * number of redshift bins
     */
    unsigned int m_Nbins;

    /**
     *  if bins are balanced or not (in terms of number of galaxies)
     */
    bool m_BalancedBins;

    /**
     * [clusters] redshift margin used to include only background galaxies
     */
    double m_ZMargin;

    /**
     * [clusters] max redshift of the halo to be kept in the cluster selection
     */
    double m_ZMaxHalo;

    /**
     * number of inpainting iterations
     */
    unsigned int m_NInpaint;

    /**
     * number of wavelet scales for inpainting
     */
    unsigned int m_NInpScale;

    /**
     * if a constraint is applied to have equal variance in/out of the masked
     * area per wavelet scale
     */
    bool m_EqualVarPerScale;

    /**
     * force B Modes to 0 in the gaps or not
     */
    bool m_ForceBMode;

    /**
     * [cartesian] if an extra border is taken into account or not
     */
    bool m_AddBorder;

    /**
     * type of projection [TAN]
     */
    std::string m_Project;

    /**
     * number of patches
     */
    unsigned int m_NPatches;

    /**
     * list of patches definitions
     */
    std::vector<PatchDef> m_Patches;

    /**
     * [sphere] HEALPix nside of the final map
     */
    int m_Nside = 0;

    /**
     * [patch2sphere] fraction of the patches overlapping
     */
    int m_overlap = 0;

    /**
     * [clusters] threshold in mass (or mass proxy: richness)
     */
    double m_MassThreshold = 0.0;

    /**
     * product extension name in fits
     */
    std::string m_ExtName;

    /**
     * parameter file type
     */
    std::string m_ParaFileType;

};  // End of GenericParam class

template<typename T>
void GenericParam::readRedshiftBaseParams(T data_xml, CatalogData& cat)
{
    // redshift parameters
    setNbins(data_xml.RedshiftBins().Nbins());
    setBalancedBins(data_xml.RedshiftBins().BalancedBins());
    auto& zbin(data_xml.RedshiftBins().RedshiftBin());
    m_ZMax.clear();
    m_ZMin.clear();
    auto catProxyNames = cat.getColumnProxyNames();
    if(not isBalancedBins())
    {
        // read bounds from params
        for (auto i = zbin.begin(); i != zbin.end(); ++i)
        {
            m_ZMax.push_back(i->ZMax());
            m_ZMin.push_back(i->ZMin());
        }
    }
    else if(std::find(catProxyNames.begin(),
                      catProxyNames.end(), "z") != catProxyNames.end() )
    {
        // catalog must be sorted by increasing redshift
        cat.sortIndex("z");
        // define redshift bin bounds manually
        double zMinParam = zbin.begin()->ZMin();
        double zMaxParam = std::prev(zbin.end())->ZMax();
        int nBins = getNbins();
        long nGal = cat.getNentriesBounds("z", zMinParam, zMaxParam);
        long nGalTot = cat.getNentries();
        double zMinBin = cat("z", 0);
        int iGal = 0;
        // if zMinParm is below lower redshift (zMinBin) start at zMinBin
        // else, search the lowest redshift
        while(zMinBin < zMinParam)
        {
            zMinBin = cat("z", iGal);
            iGal ++;
        }
        for (int i = 0; i<nBins; i++)
        {
            int nGalBin = int(nGal/nBins);
            if(nGalBin % nBins > i)
            {
                nGalBin ++;
            }
            iGal += nGalBin;
            double zMaxBin = iGal < nGalTot ? cat("z", iGal-1) : zMaxParam;
            m_ZMin.push_back(zMinBin);
            m_ZMax.push_back(zMaxBin);
            zMinBin = zMaxBin;
        }
    }
}

template<typename T>
void GenericParam::readBaseParams(T data_xml)
{
    // reduce shear parameters
    setRsNItReducedShear(data_xml.ReducedShear().NItReducedShear());
    setRsThreshold(data_xml.ReducedShear().RSThreshold());
    setRsGaussStd(data_xml.ReducedShear().GaussSTD());

    // denoising parameters
    setDenoisingAlgo(data_xml.DenoiseParams().DenoisingAlgo());
    setGaussStd(data_xml.DenoiseParams().GaussSTD());
    if (data_xml.DenoiseParams().ThresholdFDR().present())
    {
        setThresholdFdr(data_xml.DenoiseParams().ThresholdFDR().get());
    }
    else
    {
        setThresholdFdr(0);
    }

    // general parameters
    setNResamples(data_xml.NResamples());

    // data gaps parameters
    setNInpaint(data_xml.GapsParams().NInpaint());
    setEqualVarPerScale(data_xml.GapsParams().EqualVarPerScale());
    setForceBMode(data_xml.GapsParams().ForceBMode());
    setNInpScale(data_xml.GapsParams().NInpScale());
}

}  // namespace LE3_2D_MASS_WL_UTILITIES

#endif // _LE3_2D_MASS_WL_UTILITIES_GENERICPARAM_H
