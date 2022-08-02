/**
 * @file tests/src/CartesianParam_test.cpp
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

#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include "ElementsKernel/Auxiliary.h"
#include "ElementsServices/DataSync.h"

#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"

#include <ios>
#include <sstream>
#include <iostream>

using namespace LE3_2D_MASS_WL_CARTESIAN;
using namespace ElementsServices::DataSync;

namespace fs = boost::filesystem;
using boost::filesystem::path;

//-----------------------------------------------------------------------------
struct CartesianParamDataSyncFixture
{
    DataSync sync;
    path xmlFileName, clusterxmlFileName, p2sphxmlFileName, catalogFilePath,
         clusterCatalogFilePath;
    CartesianParam carParam;
    CartesianParamDataSyncFixture() :
            sync("LE3_2D_MASS_WL_CARTESIAN/datasync_webdav.conf",
                 "LE3_2D_MASS_WL_CARTESIAN/test_file_list.txt"),
            xmlFileName(sync.absolutePath("Cparam_test.xml")),
            clusterxmlFileName(sync.absolutePath("clusterParam.xml")),
            p2sphxmlFileName(sync.absolutePath("P2Sph.xml")),
            catalogFilePath(sync.absolutePath("InputLE2Catalog.xml")),
            clusterCatalogFilePath(sync.absolutePath("ClusterCatalog.xml"))
    {
        sync.download();
    }
};

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (CartesianParam_test, CartesianParamDataSyncFixture)
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE (convergencePatchxmlFile_test, CartesianParamDataSyncFixture)
{
    std::cout << "-- CARTESIANPARAM: ConvergencePatch.xml file" << std::endl;
    CatalogData cat;
    fs::path workdir = catalogFilePath.parent_path();
    fs::path catalogFileName = catalogFilePath.filename();
    cat.getCatalogData(workdir, catalogFileName);
    carParam.readConvPatchXMLFile(xmlFileName.native(), cat);
    BOOST_CHECK_EQUAL(carParam.getPatches()[0].getPatchWidth(), 10.*M_PI/180.);
    BOOST_CHECK_EQUAL(carParam.getProject(), "TAN");
    BOOST_CHECK_EQUAL(carParam.getParaFileType(), "Conv_Patch");
    BOOST_CHECK_EQUAL(carParam.getNbins(), 2);
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE (clusterxmlFile_test, CartesianParamDataSyncFixture)
{
    std::cout << "-- CARTESIANPARAM: Cluster.xml file"<< std::endl;
    CatalogData clusterCat;
    fs::path workdir = clusterCatalogFilePath.parent_path();
    fs::path catalogFileName = clusterCatalogFilePath.filename();
    clusterCat.getCatalogData(workdir, catalogFileName);
    carParam.readConvClustersXMLFile(clusterxmlFileName.native(), clusterCat);
    BOOST_CHECK_EQUAL(carParam.getPatches()[0].getPatchWidth(), 10.*M_PI/180.);
    BOOST_CHECK_EQUAL(carParam.getZMaxHalo(), 0.6);
    BOOST_CHECK_EQUAL(carParam.getZMargin(), 0.0);
    BOOST_CHECK_EQUAL(carParam.getMassThreshold(), 0.1);
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE (patches2SpherexmlFile_test, CartesianParamDataSyncFixture)
{
    std::cout << "-- CARTESIANPARAM: Patches2Sphere.xml file"<< std::endl;
    double placeHolder = 0.;
    CatalogData dummy;
    carParam.readConvPatchesToSphereXMLFile(p2sphxmlFileName.native(),
                placeHolder, placeHolder, placeHolder, placeHolder, dummy);
    BOOST_CHECK_EQUAL(carParam.getNside(), 512);
}

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(readParameterFile_test, CartesianParamDataSyncFixture)
{
    std::cout << "-- CARTESIANPARAM: readParameterFile_test"<< std::endl;
    CatalogData dummy;
    readParameterFile(xmlFileName, carParam, dummy);
    BOOST_CHECK_EQUAL(carParam.getPatches()[0].getPatchWidth(), 10.*M_PI/180.);
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(readParameterFileCluster_test, CartesianParamDataSyncFixture)
{
    std::cout << "-- CARTESIANPARAM: readParameterFileCluster_test"<< std::endl;
    CatalogData dummy;
    CatalogData clusterCat;
    fs::path workdir = clusterCatalogFilePath.parent_path();
    fs::path catalogFileName = clusterCatalogFilePath.filename();
    clusterCat.getCatalogData(workdir, catalogFileName);
    readParameterFile(clusterxmlFileName, carParam, dummy, 0,0,0,0, clusterCat);
    BOOST_CHECK_EQUAL(carParam.getPatches()[0].getPatchWidth(), 10.*M_PI/180.);
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(readParameterFilePatches2Sphere_test, CartesianParamDataSyncFixture)
{
    std::cout << "-- CARTESIANPARAM: readParameterFilePatches2Sphere_test"<< std::endl;
    CatalogData dummy;
    readParameterFile(p2sphxmlFileName, carParam, dummy, 70., 90., 350., 360.);
    BOOST_CHECK_EQUAL(carParam.getNside(), 512);
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
