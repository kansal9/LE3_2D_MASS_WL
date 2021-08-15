/**
 * @file src/lib/SplitCatalog.cpp
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

#include "LE3_2D_MASS_WL_CARTESIAN/SplitCatalog.h"

using namespace Euclid::WeakLensing::TwoDMass;

static Elements::Logging logger = Elements::Logging::getLogger("Splitting");
namespace LE3_2D_MASS_WL_CARTESIAN {

SplitCatalog::SplitCatalog(std::vector<std::vector<double> >& Cat_Data, const CartesianParam& param,
          std::string& catType): InData(Cat_Data), m_cartesianParam(param), m_catType(catType) {
   m_nbZBins = m_cartesianParam.getnbZBins();
   m_balancedBins = m_cartesianParam.get_BalancedBins();
   double redshiftMin, redshiftMax;
   vecMinMax(InData[5], &redshiftMin, &redshiftMax);
   std::vector<double> zMin = m_cartesianParam.getZMin();
   double zMax = m_cartesianParam.getZMax();
   if (zMin[0] < redshiftMin) {
    m_zMin = redshiftMin;
   } else {
    m_zMin = zMin[0];
   }
   if (zMax > redshiftMax) {
    m_zMax = redshiftMax;
   } else {
    m_zMax = zMax;
   }
}

SplitCatalog::SplitCatalog(std::vector<std::vector<double> >& Cat_Data,
                    const LE3_2D_MASS_WL_SPHERICAL::SphericalParam& Sphparam, std::string& catType):
                                       InData(Cat_Data), m_SphericalParam(Sphparam), m_catType(catType) {
   m_nbZBins = m_SphericalParam.getNbins();
   m_balancedBins = m_SphericalParam.getBalancedBins();
   double redshiftMin, redshiftMax;
   vecMinMax(InData[5], &redshiftMin, &redshiftMax);
   if (m_SphericalParam.getZMin() < redshiftMin) {
    m_zMin = redshiftMin;
   } else {
    m_zMin = m_SphericalParam.getZMin();
   }
   if (m_SphericalParam.getZMax() > redshiftMax) {
    m_zMax = redshiftMax;
   } else {
    m_zMax = m_SphericalParam.getZMax();
   }
}
/*
void SplitCatalog::getColname (std::vector < std::string >& columnNames) {

    const std::string ra = "SHE_" + m_catType + "_UPDATED_RA";
    const std::string dec = "SHE_" + m_catType + "_UPDATED_DEC";
    const std::string g1= "SHE_" + m_catType + "_G1";
    const std::string g2 = "SHE_" + m_catType + "_G2";
    const std::string z = "PHZ_MEDIAN";
    const std::string w = "SHE_" + m_catType + "_WEIGHT";
    const std::string phzCorrection = "PHZ_" + m_catType + "_CORRECTION";
    columnNames.push_back(ra);
    columnNames.push_back(dec);
    columnNames.push_back(g1);
    columnNames.push_back(g2);
    columnNames.push_back(z);
    columnNames.push_back(w);
    columnNames.push_back(phzCorrection);
}*/

void SplitCatalog::getEqualBinSplittedCatalogs(fs::path& datadir, std::vector<std::string>& Filenames) {

   std::vector < std::string > columnNames = getcolumnNames(m_catType);
   //getColname (columnNames);

   for (size_t iter = 0; iter<m_nbZBins; iter++) {
    const std::string filename = datadir.string() + "/" + Filenames[iter];
    std::vector <double> temp;
    temp = split(InData[0], m_nbZBins, iter);
    saveSubCatalogs (filename, temp, columnNames[0]);
    temp = split(InData[1], m_nbZBins, iter);
    saveSubCatalogs (filename, temp, columnNames[1]);
    temp = split(InData[3], m_nbZBins,  iter);
    saveSubCatalogs (filename, temp, columnNames[2]);
    temp = split(InData[4], m_nbZBins, iter);
    saveSubCatalogs (filename, temp, columnNames[3]);
    temp = split(InData[5], m_nbZBins, iter);
    saveSubCatalogs (filename, temp, columnNames[4]);
    temp = split(InData[6], m_nbZBins, iter);
    saveSubCatalogs (filename, temp, columnNames[5]);
    temp = split(InData[7], m_nbZBins, iter);
    saveSubCatalogs (filename, temp, columnNames[6]);
  }
}

void SplitCatalog::getSplittedCatalogs(fs::path& datadir, std::vector<std::string>& Filenames) {
  // unsigned int nGal=0;
  // nGal = InData[0].size();
   //logger.info()<<"number of galaxies in input Catalog: "<< nGal;
   double range = (m_zMax-m_zMin)/m_nbZBins;
   double t_zmin = m_zMin;
   double t_zmax= m_zMin + range;
   //logger.info()<<"zMax: "<< m_zMax <<"\t"<<"zMin: "<<m_zMin<<"\t"<<"range: "<<range;
   std::vector < std::string > columnNames = getcolumnNames(m_catType);
   //getColname (columnNames);
   for (size_t it = 0; it<m_nbZBins; it++) {
    const std::string filename = datadir.string() + "/" + Filenames[it];
    std::vector <int> indices;
    indices.clear();
    logger.info()<<"Redshift Range for subcatalogue " << int(it+1) <<"is z_min: "<< t_zmin <<"\t"<<"z_max: "<<t_zmax;
    if (t_zmax > m_zMax) {
     logger.info() <<" Redshift values are greater than input Zmax ";
     break;
    }
 /*   for (size_t i=0; i<nGal; i++){
      if (t_zmax<=m_zMax) {
        if (InData[5][i] >= t_zmin && InData[5][i] < t_zmax) { //less than max and greater than and equal to min
         indices.push_back(i);
        }
      }
      if ((m_nbZBins % 2 == 1) && t_zmax > m_zMax) { //left overs
        if (InData[5][i] >= t_zmin && InData[5][i] <= m_zMax){
         indices.push_back(i);
        }
      }
   }*/
   getIndices(indices, t_zmin, t_zmax);
   indices.shrink_to_fit();
   std::vector <double> temp;
   temp = extract_at(InData[0], indices);
   saveSubCatalogs (filename, temp, columnNames[0]);
   temp = extract_at(InData[1], indices);
   saveSubCatalogs (filename, temp, columnNames[1]);
   temp = extract_at(InData[3], indices);
   saveSubCatalogs (filename, temp, columnNames[2]);
   temp = extract_at(InData[4], indices);
   saveSubCatalogs (filename, temp, columnNames[3]);
   temp =  extract_at(InData[5], indices);
   saveSubCatalogs (filename, temp, columnNames[4]);
   temp = extract_at(InData[6], indices);
   saveSubCatalogs (filename, temp, columnNames[5]);
   temp = extract_at(InData[7], indices);
   saveSubCatalogs (filename, temp, columnNames[6]);

   logger.info()<<"number of created subcatalogs: "<<int(it+1);
   t_zmin = t_zmin + range;
   t_zmax = t_zmax + range;
  }
}

void SplitCatalog::getIndices(std::vector < int >& Indices, double t_zmin, double t_zmax) {
      unsigned int nGal=0;
      nGal = InData[0].size();
      for (size_t i=0; i<nGal; i++){
      if (t_zmax<=m_zMax) {
        if (InData[5][i] >= t_zmin && InData[5][i] < t_zmax) { //less than max and greater than and equal to min
         Indices.push_back(i);
        }
      }
      if ((m_nbZBins % 2 == 1) && t_zmax > m_zMax) { //left overs
        if (InData[5][i] >= t_zmin && InData[5][i] <= m_zMax){
         Indices.push_back(i);
        }
      }
   }
}

bool SplitCatalog::writeSubCatalogs(fs::path& datadir, std::vector<std::string>& Filenames){
 logger.info()<<"number of redshift bins are:"<< m_nbZBins;
 logger.info()<<"Balanced bins:"<< m_balancedBins;

  if (m_zMin < 0.0 ){
   logger.info()<< "Warning: Redshift Value is less than ZERO. The minimum redshift value is: "<< m_zMin;
  }
  if (m_zMax > 10.0){
   logger.info()<< "Warning: Redshift Value is greater than TEN. The maximum redshift value is: "<< m_zMax;
  }
   if (InData[2].empty()) {
    InData[2].resize(InData[5].size());
   }

  if (1 == m_balancedBins) {
   logger.info()<<"entering Balanced bins";
   getEqualBinSplittedCatalogs(datadir, Filenames);
  } else { //end of loop for balanced bins and start for Zbin split
   logger.info()<<"entering split";
   getSplittedCatalogs(datadir, Filenames);
  }

return true;
}

bool SplitCatalog::saveSubCatalogs (const std::string& filename, std::vector<double>& Data, const std::string &colname){
  if (access( filename.c_str(), F_OK ) != -1 ) {
    MefFile outfile(filename, MefFile::Permission::Edit);
    const auto col = generateColumn<double>(colname, Data);
    long nbHdu = outfile.hduCount();
    const auto &ext = outfile.access<BintableHdu>(nbHdu-1);
    ext.appendColumn(col);
  } else {
    MefFile outfile(filename, MefFile::Permission::Overwrite);
    const auto col = generateColumn<double>(colname, Data);
    outfile.assignBintableExt("", col); // Unnamed extension
  }

 return true;
}

}  // namespace LE3_2D_MASS_WL_CARTESIAN
