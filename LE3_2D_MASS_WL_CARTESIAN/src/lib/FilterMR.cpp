/**
 * @file src/lib/FilterMR.cpp
 * @date 03/12/21
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

#include "LE3_2D_MASS_WL_CARTESIAN/FilterMR.h"

using namespace LE3_2D_MASS_WL_UTILITIES;

static Elements::Logging logger = Elements::Logging::getLogger("FilterMR");

namespace LE3_2D_MASS_WL_CARTESIAN
{

FilterMR::FilterMR(ConvergenceMap &convMap, bool PositiveCons,
        bool KillLastScale, bool KillIsol, int niter, int FirstScale,
        LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam) :
        m_convMap(convMap), m_KillLastScale(KillLastScale), m_KillIsol(
                KillIsol), m_cartesianParam(cartesianParam), m_MP(
                convMap.getXdim(), convMap.getYdim())

{

    // Retrieve the size of the maps
    m_Xaxis = convMap.getXdim();
    m_Yaxis = convMap.getYdim();
    m_Zaxis = convMap.getZdim();
    m_PositveConstraint = PositiveCons;
    m_FirstScale = FirstScale;
    m_niter = niter;
    m_MaxIter = 20;
    m_SigmaNoise = 0.0;
    m_avar = 2.0;
    logger.info() << "axis dim: " << m_Xaxis << " " << m_Yaxis << " "
            << m_Zaxis;
    logger.info() << std::boolalpha << "Kill Last Scale = " << m_KillLastScale;
    logger.info() << std::boolalpha << "kill isolated pixels = " << m_KillIsol;
    logger.info() << std::boolalpha << "apply positive constraint = "
            << m_PositveConstraint;
    logger.info() << "First Scale = " << m_FirstScale - 1;
    m_FDR = m_cartesianParam.getThresholdFdr();
    if (m_FDR == 0)
    {
        m_FDR = 0.05;
    }
    logger.info() << "FDR_val = " << m_FDR;
    m_nbScales = m_cartesianParam.getNInpScale();
    if (m_nbScales == 0)
    {
        m_nbScales = int(log(m_Xaxis) / log(2.)) - 2.;
    }
    logger.info() << "nbScales: " << m_nbScales;
}

FilterMR::FilterMR(const FilterMR& copy) :
        m_convMap(copy.m_convMap), m_Xaxis(copy.m_Xaxis), m_Yaxis(copy.m_Yaxis), m_Zaxis(
                copy.m_Zaxis), m_nbScales(copy.m_nbScales), m_niter(
                copy.m_niter), m_FirstScale(copy.m_FirstScale), m_MaxIter(
                copy.m_MaxIter), m_SigmaNoise(copy.m_SigmaNoise), m_FDR(
                copy.m_FDR), m_avar(copy.m_avar), m_PositveConstraint(
                copy.m_PositveConstraint), m_KillLastScale(
                copy.m_KillLastScale), m_KillIsol(copy.m_KillIsol), m_cartesianParam(
                copy.m_cartesianParam), m_MP(copy.m_MP)
{
}

double FilterMR::fdr_pv(std::vector<double> &PValue, double alpha)
{
    std::vector<double> input_pval;
    for (size_t i = 0; i < PValue.size(); i++)
    {
        input_pval.push_back(PValue[i]);
    }
    std::sort(input_pval.begin(), input_pval.end());
//    input_pval[1] =0;
    int m = input_pval.size();
    double p_cutoff = 0.;
    for (int i = 0; i < m; i++)
    {
        double a = alpha * i / m;
        if (input_pval[i] < a)
        {
            p_cutoff = input_pval[i];
        }
    }
    return p_cutoff;
}

double FilterMR::myErfInv(double x)
{
    double tt1, tt2, lnx, sgn;
    sgn = (x < 0) ? -1.0 : 1.0;
    x = (1 - x) * (1 + x); // x = 1 - x*x;
    lnx = log(x);

    tt1 = 2 / (M_PI * 0.147) + 0.5 * lnx;
    tt2 = 1 / (0.147) * lnx;

    return (sgn * sqrt(-tt1 + sqrt(tt1 * tt1 - tt2)));
}

double FilterMR::xerfc(double F)
{
    double Nu = 0;
    double P = F;
    if (P >= 1.)
    {
        P = 1.;
    }
    else
    {
        if (P < 0)
        {
            P = 0.;
        }
        if (P > 0.5)
        {
            Nu = sqrt(2.) * myErfInv(1 - ((1. - P) / 0.5));
        }
        else
        {
            Nu = -sqrt(2.) * xerfc((double) P / 0.5);
        }
    }
    return Nu;
}

double FilterMR::getProbility(double val, double sigLevel)
{
    double P = 0., vn = 0.;
    if (fabs(val) < FLOAT_EPSILON)
    {
        P = 1.;
    }
    else
    {
        if (sigLevel < FLOAT_EPSILON)
        {
            P = 0.;
        }
        else
        {
            vn = fabs(val) / (sqrt(2.) * sigLevel);
            if (vn > 3.5)
            {
                P = 0.;
            }
            else
            {
                P = std::erfc(vn);
            }
        }
    }
    return P;
}

void FilterMR::get_fdr(std::vector<Matrix>& myBand, std::vector<double>& NSigma)
{
    logger.info() << "Begin";

    for (int kScale = 0; kScale < m_nbScales - 1; kScale++)
    {
        std::vector<double> PValue;
        for (int j = 0; j < m_Yaxis; j++)
        {
            for (int i = 0; i < m_Xaxis; i++)
            {
                PValue.push_back(
                        getProbility(myBand[kScale](i, j),
                                (m_SigmaNoise * NormB3Spline[kScale])));
            }
        }
        double alpha_b = m_FDR * pow(m_avar, double(kScale));
        //alpha_b = (kScale==0) ? 1. - erf(NSigma[kScale] / sqrt(2.)): alpha_b*2;
        if (alpha_b > 0.5)
        {
            alpha_b = 0.5;
        }
        double Pdet = fdr_pv(PValue, alpha_b);

        double temp = 0.5 + (1. - Pdet) / 2.;
        NSigma[kScale] = fabs(xerfc(temp));
        if ((NSigma[kScale] < 7) && (NSigma[kScale] > 0))
        {
            NSigma[kScale] = fabs(xerfc(temp));
        }
        else
        {
            NSigma[kScale] = 7.;
        }

        logger.info() << "band: " << kScale + 1 << "    m_FDR: " << alpha_b
                << "    nsigma: " << NSigma[kScale];
    }

    logger.info() << "End";
}

void FilterMR::kill_isolated(Matrix& image)
{
    Matrix has_neighbor(m_Xaxis, m_Yaxis);
    has_neighbor.clear(); // set all values to 0

    for (int j = 1; j < m_Yaxis - 1; j++)
    {
        for (int i = 1; i < m_Xaxis - 1; i++)
        {
            if ((has_neighbor(i, j) || image(i - 1, j) > 0)
                    && (has_neighbor(i, j) || image(i + 1, j) > 0)
                    && (has_neighbor(i, j) || image(i, j - 1) > 0)
                    && (has_neighbor(i, j) || image(i, j + 1) > 0))
            {
                has_neighbor(i, j) = 1;
            }
        }
    }
    for (int j = 0; j < m_Yaxis; j++)
    {
        for (int i = 0; i < m_Xaxis; i++)
        {
            if (has_neighbor(i, j) == 1)
            {
                image(i, j) = 0;
            }
        }
    }
}

void FilterMR::set_support(std::vector<Matrix>& myBand,
        std::vector<double>& NSigma, int it)
{
    for (int j = 0; j < m_Yaxis; j++)
    {
        for (int i = 0; i < m_Xaxis; i++)
        {
            double temp = myBand[it](i, j);
            if (temp
                    < (NSigma[it] * m_SigmaNoise * NormB3Spline[it]))
            {
                myBand[it](i, j) = 0;
            }
            if (temp > 0)
            {
                myBand[it](i, j) = 1;
            }
        }
    }
}

void FilterMR::mr_support(std::vector<Matrix>& myBand,
        std::vector<double>& NSigma, std::vector<Matrix>& mr_myBand)
{

    std::vector<Matrix> copy_myBand = myBand;
    for (int it = m_FirstScale - 1; it < m_nbScales - 1; it++)
    {
        set_support(copy_myBand, NSigma, it);

        if (m_KillIsol == true)
        {
            kill_isolated(copy_myBand[it]);
        }
        for (int j = 0; j < m_Yaxis; j++)
        {
            for (int i = 0; i < m_Xaxis; i++)
            {
                mr_myBand[it](i, j) = copy_myBand[it](i, j);
            }
        }
    }
}

void FilterMR::performFiltering(ConvergenceMap& outConvergenceMap)
{
    // Will work only on kappaE (first axis k = 0 of the GenericMap)

    // Create a copy of the conv map
    ConvergenceMap kappaMap(m_convMap);

    m_SigmaNoise = kappaMap.getSigma(0);
    logger.info() << "sigmaNoise: " << m_SigmaNoise;

    std::vector<Matrix> myBand;
    m_MP.transformBspline(kappaMap[0], myBand, m_nbScales);

    std::vector<double> NSigma_l(m_nbScales - 1);

    get_fdr(myBand, NSigma_l);

    std::vector<Matrix> mr_myBand = myBand;
    for (size_t b = 0; b < mr_myBand.size(); b++)
    {
        mr_myBand[b].clear();
    }

    mr_support(myBand, NSigma_l, mr_myBand);

    mr_myBand[m_nbScales - 1].clear();

    applyfilter(myBand, NSigma_l, outConvergenceMap);

    outConvergenceMap.clear(0);

    reconstruct_optimized(myBand, mr_myBand, outConvergenceMap);

    if (true == m_PositveConstraint)
    {
        outConvergenceMap.applyThreshold(0.);
    }
}

void FilterMR::reconstruct_optimized(std::vector<Matrix>& band,
        std::vector<Matrix>& mr_band, ConvergenceMap& outConvergenceMap)
{

    for (int b = 0; b < m_FirstScale - 1; b++)
    {
        band[b].clear();
    }

    //Kill last scale
    if (true == m_KillLastScale)
    {
        band[m_nbScales - 1].clear();
    }

    reconstruct(band, outConvergenceMap[0]);

    std::vector<Matrix> band_i;
    std::vector<Matrix> band_res;

    for (int it = 0; it < m_niter; it++)
    {
        std::cout << "start iteration " << it + 1 << std::endl;
        m_MP.transformBspline(outConvergenceMap[0], band_i, m_nbScales);
        band_res = band_i;
        for (int b = 0; b < m_nbScales; b++)
        {
            for (int j = 0; j < m_Yaxis; j++)
            {
                for (int i = 0; i < m_Xaxis; i++)
                {
                    double a0 = band[b](i, j);
                    double b0 = band_i[b](i, j);
                    band_res[b](i, j) = a0 - b0;
                    double c0 = mr_band[b](i, j);
                    if (c0 == 0)
                    {
                        band_res[b](i, j) = 0.;
                    }
                }
            }
        }

        for (int j = 0; j < m_Yaxis; j++)
        {
            for (int i = 0; i < m_Xaxis; i++)
            {
                double a1 = band[m_nbScales - 1](i, j);
                double b1 = band_i[m_nbScales - 1](i, j);
                band_res[m_nbScales - 1](i, j) = a1 - b1;
            }
        }

        Matrix res(m_Xaxis, m_Yaxis);
        reconstruct(band_res, res);

        for (int j = 0; j < m_Yaxis; j++)
        {
            for (int i = 0; i < m_Xaxis; i++)
            {
                double a2 = outConvergenceMap.getBinValue(i, j, 0);
                double b2 = res(i, j);
                if (b2 > 0)
                {
                    outConvergenceMap.setBinValue(i, j, 0, a2 + b2);
                }
            }
        }
    }
}

void FilterMR::reconstruct(std::vector<Matrix>& band, Matrix& outConvergenceMap)
{
    for (int j = 0; j < m_Yaxis; j++)
    {
        for (int i = 0; i < m_Xaxis; i++)
        {
            outConvergenceMap(i, j) = band[m_nbScales - 1](i, j);
        }
    }

    Matrix temp(m_Xaxis, m_Yaxis);

    size_t s = m_nbScales - 2;
    do
    {
        m_MP.smoothBspline(outConvergenceMap, temp, s);
        for (int j = 0; j < m_Yaxis; j++)
        {
            for (int i = 0; i < m_Xaxis; i++)
            {
                double aa = temp(i, j);
                double bb = band[s](i, j);
                outConvergenceMap(i, j) = aa + bb;
            }
        }

    } while (s--);
}

void FilterMR::applyfilter(std::vector<Matrix>& band,
        std::vector<double>& NSigma, ConvergenceMap& outConvergenceMap)
{
    std::vector<double> TabAlpha(m_nbScales - 1);
    std::fill(TabAlpha.begin(), TabAlpha.end(), 1.);
    if(TabAlpha.size() > 0)
    {
        TabAlpha[0] = 5.;
    }
    std::vector<double> TabDelta(m_nbScales - 1);
    std::vector<double> regulMin(m_nbScales - 1);
    std::fill(regulMin.begin(), regulMin.end(), 0.);
    std::vector<double> regulMax(m_nbScales - 1);
    std::fill(regulMax.begin(), regulMax.end(), 100.);

    int iter = 0;

    std::vector<Matrix> copyBand;

    ReconstructMR restore;

    do
    {
        iter++;
        logger.info() << " START Multiscale Entropy Iteration = " << iter;

        copyBand = band;

        for (int b = 0; b < m_nbScales - 1; b++)
        {
            TabAlpha[b] = (regulMin[b] + regulMax[b]) / 2.;
        }

        fixedAlphaFilter(copyBand, TabAlpha, NSigma, restore);

        m_MP.reconsBspline(copyBand, outConvergenceMap[0]);

        double Flux = outConvergenceMap.getFlux();
        double Mean = Flux / (m_Yaxis * m_Xaxis);

        logger.info() << " Min = " << outConvergenceMap.getMin() << "    Max = "
                << outConvergenceMap.getMax() << "    Flux = " << Flux;
        logger.info() << "    Mean = " << Mean << "     Sigma = "
                << outConvergenceMap.getSigma();

        if (true == m_PositveConstraint)
        {
            outConvergenceMap.applyThreshold(0.);
        }

        m_MP.transformBspline(outConvergenceMap[0], copyBand, m_nbScales);
        estimateNewAlphaParam(band, copyBand, regulMin, regulMax, TabAlpha,
                TabDelta);
        logger.info() << " END Multiscale Entropy Iteration = " << iter;

    } while (((fabs(*max_element(TabDelta.begin(), TabDelta.end()))) > 0.01)
            && (iter < m_MaxIter));

    for (int b = 0; b < m_nbScales - 1; b++)
    {
        TabAlpha[b] = (regulMin[b] + regulMax[b]) / 2;
    }

    fixedAlphaFilter(band, TabAlpha, NSigma, restore);

}

void FilterMR::estimateNewAlphaParam(std::vector<Matrix>& band,
        std::vector<Matrix>& band_sol, std::vector<double>& regulMin,
        std::vector<double>& regulMax, std::vector<double>& TabAlpha,
        std::vector<double>& TabDelta)
{
    double sigmaNoise;
    for (int b = m_FirstScale - 1; b < m_nbScales - 1; b++)
    {
        sigmaNoise = 0.;
        double sigmaCoeff = m_SigmaNoise * NormB3Spline[b];

        for (int j = 0; j < m_Yaxis; j++)
        {
            for (int i = 0; i < m_Xaxis; i++)
            {
                double residual = band[b](j, i) - band_sol[b](j, i);
                double var = sigmaCoeff * sigmaCoeff;
                sigmaNoise = sigmaNoise + ((residual * residual) / var);
            }
        }
        sigmaNoise = sqrt(sigmaNoise / (double) (m_Yaxis * m_Xaxis));
        if (sigmaNoise >= 1)
        {
            regulMax[b] = TabAlpha[b];
        }
        else
        {
            regulMin[b] = TabAlpha[b];
        }
        TabDelta[b] = fabs(sigmaNoise - 1.);
    }
}

void FilterMR::fixedAlphaFilter(std::vector<Matrix>& band,
        std::vector<double>& TabAlpha, std::vector<double>& NSigma,
        LE3_2D_MASS_WL_CARTESIAN::ReconstructMR & restore)
{

    double alphaP, val = 0.;
    for (int b = 0; b < m_nbScales - 1; b++)
    {
        //  std::cout << "rp: " << TabAlpha[b] <<std::endl;
        double sigma = m_SigmaNoise * NormB3Spline[b];
        for (int j = 0; j < m_Yaxis; j++)
        {
            for (int i = 0; i < m_Xaxis; i++)
            {
                if (b < m_FirstScale - 1)
                {
                    band[b](i, j) = 0.;
                }
                else
                {
                    double rp = TabAlpha[b];
                    val = band[b](i, j);
                    if (val >= 0)
                    {
                        alphaP = fabs(val / (sigma * NSigma[b]));
                        if (alphaP > 1)
                        {
                            alphaP = 1;
                        }
                        if (alphaP < FLOAT_EPSILON)
                        {
                            rp = -1;
                        }
                        if (alphaP != 0)
                        {
                            rp *= (1 - alphaP) / alphaP;
                        }
                    }
                    // std::cout << "regularParam: " << rp << std::endl;
                    double temp_val = restore.filter(val, rp, sigma);
                    if (rp > 0)
                    {
                        band[b](i, j) = temp_val;
                    }
                    else
                    {
                        band[b](i, j) = val;
                    }
                }
            }
        }
    }
}

}  // namespace LE3_2D_MASS_WL_CARTESIAN
