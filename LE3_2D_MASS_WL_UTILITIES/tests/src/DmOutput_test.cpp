/**
 * @file tests/src/DmOutput_test.cpp
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

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmOutput.h"

#include "ElementsKernel/Temporary.h"
#include <fstream>
#include <string>

using namespace dpd::le3::wl::twodmass::out::convergencepatch;
using namespace dpd::le3::wl::twodmass::out::convergencesphere;
using namespace dpd::le3::wl::twodmass::out::convergencepatchestosphere;
using namespace pro::le3::wl::twodmass;
using namespace LE3_2D_MASS_WL_UTILITIES;

namespace fs = boost::filesystem;

using Elements::TempDir;
using Elements::TempFile;

using std::string;

//----------------------------------------------------------------------------
struct DmOutputTestFixture
{
    DmOutput dm;
    TempDir td;
    fs::path test_path;

    DmOutputTestFixture() : dm(), td()
    {
        test_path = td.path();
        // test_path = "/tmp";
        std::cout << "test_path: " << test_path << std::endl;
    }
};


//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE (DmOutput_test)
//-----------------------------------------------------------------------------
// BOOST_FIXTURE_TEST_CASE(getGenericHeader_test, DmOutputTestFixture)
BOOST_AUTO_TEST_CASE(getGenericHeader_test)
{
    try {
        // Create generator
        const std::string productType = "TestProduct";
        GenericHeaderGenerator generator(productType);
        sys::genericHeader* header = generator.generate();
        delete header;
    } catch (std::exception& e) {
        BOOST_THROW_EXCEPTION(e);
    }
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(createPatchXml_test, DmOutputTestFixture)
{
    std::cout << "-- DmOutput:createPatchXml_test" << std::endl;
    try
    {
        std::string product_type = "DpdTwoDMassConvergencePatch";
        int NResamples = 0;
        auto product = initProduct<dpdTwoDMassConvergencePatch,
                                   twoDMassCollectConvergencePatch,
                                   int>(product_type, NResamples);

        TempFile f1, f2, f3, f4;
        fs::ofstream ofs(f1.path());
        ofs.close();
        ofs.open(f2.path());
        ofs.close();
        ofs.open(f3.path());
        ofs.close();

        const fs::path file = test_path / "Patch.xml";

        // add existing file
        dm.createPatchXml(product, outputType::NoisedPatch, f1.path());
        dm.createPatchXml(product, outputType::DenoisedPatch, f2.path());
        dm.createPatchXml(product, outputType::SNRPatch, f3.path());

        // nothing will be done if file does not exist
        dm.createPatchXml(product, outputType::SNRPatch, f4.path());

        writeProduct<dpdTwoDMassConvergencePatch>(product, file);

    BOOST_CHECK(fs::is_regular_file(file));
    } catch (std::exception& e) {
        BOOST_THROW_EXCEPTION(e);
    }
}
//-----------------------------------------------------------------------------
/*
BOOST_FIXTURE_TEST_CASE(createSingleClusterXml_test, DmOutputTestFixture)
{
    std::cout << "-- DmOutput:createSingleClusterXml_test" << std::endl;
    try
    {
        fs::path fitsFile{ "testFile.fits" };
        const fs::path filename = test_path / "ClusterConvergenceMapPatch.xml";
        std::string product_type = "DpdTwoDMassConvergenceClusters";
        std::string ClusterID = "198";

        auto product = initProduct<dpdTwoDMassConvergenceSingleCluster,
                                   twoDMassConvergenceSingleCluster,
                                   std::string>(product_type, ClusterID);

        dm.createSingleClusterXml(product, fitsFile);
        writeProduct<dpdTwoDMassConvergenceSingleCluster>(product, filename);
        BOOST_CHECK(fs::is_regular_file(filename));
    } catch (std::exception& e) {
        BOOST_THROW_EXCEPTION(e);
    }
}
*/
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(createSphereOutputXml_test, DmOutputTestFixture)
{
    std::cout << "-- DmOutput:createSphereOutputXml_test" << std::endl;
    try
    {
        fs::path out_fits_file_1{ "some_file_1.fits" };
        fs::path out_fits_file_2{ "some_file_2.fits" };
        fs::path out_fits_file_3{ "some_file_3.fits" };
        fs::path out_fits_file_4{ "some_file_4.fits" };
        fs::path out_fits_file_5{ "some_file_5.fits" };
        fs::path out_fits_file_6{ "some_file_6.fits" };
        const fs::path filename = test_path / "Sphere.xml";
        std::string product_type = "DpdTwoDMassConvergenceSphere";
        int NResamples = 2;

        auto product = initProduct<dpdTwoDMassConvergenceSphere,
                                   twoDMassCollectConvergenceSphere,
                                   int>(product_type, NResamples);

        dm.createNoisedSphereXml(product, out_fits_file_1);
        dm.createDenoisedSphereXml(product, out_fits_file_2);
        dm.createSNRSphereOutputXml(product, out_fits_file_3);
        dm.createSphereGalCountXml(product, out_fits_file_4);
        dm.createSphereMCXml(product, out_fits_file_5);
        dm.createSphereMCXml(product, out_fits_file_6);

        writeProduct<dpdTwoDMassConvergenceSphere>(product, filename);
        BOOST_CHECK(fs::is_regular_file(filename));
    } catch (std::exception& e) {
        BOOST_THROW_EXCEPTION(e);
    }
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(createNoisedSphereXml_test, DmOutputTestFixture)
{
    std::cout << "-- DmOutput:createNoisedSphereXml_test" << std::endl;
    try
    {
        fs::path fitsFileSphere{ "testFile.fits" };
        const fs::path filename = test_path / "NoisedConvergenceSphere.xml";
        std::string product_type = "DpdTwoDMassConvergenceSphere";
        int NResamples = 0;

        auto product =  initProduct<dpdTwoDMassConvergenceSphere,
                                    twoDMassCollectConvergenceSphere,
                                    int>(product_type, NResamples);

        dm.createNoisedSphereXml(product, fitsFileSphere);
        writeProduct<dpdTwoDMassConvergenceSphere>(product, filename);
        BOOST_CHECK(fs::is_regular_file(filename));
    } catch (std::exception& e) {
        BOOST_THROW_EXCEPTION(e);
    }
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(createDenoisedSphereXml_test, DmOutputTestFixture)
{
    std::cout << "-- DmOutput:createDenoisedSphereXml_test" << std::endl;
    try
    {
        fs::path fitsFileSphere{ "testFile.fits" };
        const fs::path filename = test_path / "DeNoisedConvergenceSphere.xml";
        std::string product_type = "DpdTwoDMassConvergenceSphere";
        int NResamples = 0;

        auto product = initProduct<dpdTwoDMassConvergenceSphere,
                                   twoDMassCollectConvergenceSphere,
                                   int>(product_type, NResamples);

        dm.createDenoisedSphereXml(product, fitsFileSphere);
        writeProduct<dpdTwoDMassConvergenceSphere>(product, filename);
        BOOST_CHECK(fs::is_regular_file(filename));
    } catch (std::exception& e) {
        BOOST_THROW_EXCEPTION(e);
    }
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(createSphereGalCountXml_test, DmOutputTestFixture)
{
    std::cout << "-- DmOutput:createSphereGalCountXml_test" << std::endl;
    try
    {
        fs::path fitsFileSphere{ "testFile.fits" };
        const fs::path filename = test_path / "GalCountSphere.xml";
        std::string product_type = "DpdTwoDMassConvergenceSphere";
        int NResamples = 0;

        auto product = initProduct<dpdTwoDMassConvergenceSphere,
                                   twoDMassCollectConvergenceSphere,
                                   int>(product_type, NResamples);

        dm.createSphereGalCountXml(product, fitsFileSphere);
        writeProduct<dpdTwoDMassConvergenceSphere>(product, filename);
        BOOST_CHECK(fs::is_regular_file(filename));
    } catch (std::exception& e) {
        BOOST_THROW_EXCEPTION(e);
    }
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(createSphereMCXml_test, DmOutputTestFixture)
{
    std::cout << "-- DmOutput:createSphereMCXml_test" << std::endl;
    try
    {
        fs::path fitsFileSphere{ "testFile.fits" };
        const fs::path filename = test_path / "MCConvergenceSphere.xml";
        std::string product_type = "DpdTwoDMassConvergenceSphere";
        int NResamples = 1;

        auto product = initProduct<dpdTwoDMassConvergenceSphere,
                                   twoDMassCollectConvergenceSphere,
                                   int>(product_type, NResamples);

        dm.createSphereMCXml(product, fitsFileSphere);
        writeProduct<dpdTwoDMassConvergenceSphere>(product, filename);
        BOOST_CHECK(fs::is_regular_file(filename));
    } catch (std::exception& e) {
        BOOST_THROW_EXCEPTION(e);
    }
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(createSNRSphereOutputXml_test, DmOutputTestFixture)
{
    std::cout << "-- DmOutput:createSNRSphereOutputXml_test" << std::endl;
    try
    {
        fs::path fitsFileSphere{ "testFile.fits" };
        const fs::path filename = test_path / "SNRSphere.xml";
        std::string product_type = "DpdTwoDMassConvergenceSphere";
        int NResamples = 2;

        auto product = initProduct<dpdTwoDMassConvergenceSphere,
                                   twoDMassCollectConvergenceSphere,
                                   int>(product_type, NResamples);

        dm.createSNRSphereOutputXml(product, fitsFileSphere);
        writeProduct<dpdTwoDMassConvergenceSphere>(product, filename);
        BOOST_CHECK(fs::is_regular_file(filename));
    } catch (std::exception& e) {
        BOOST_THROW_EXCEPTION(e);
    }
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(createPeakOutputXml_test, DmOutputTestFixture)
{
    std::cout << "-- DmOutput: createPeakOutputXml_test" << std::endl;
    try
    {
        fs::path fitsFilePC{ "testFile.fits" };
        const fs::path filename = test_path / "PeakCount.xml";
        dm.createPeakOutputXml(filename, fitsFilePC);
        BOOST_CHECK(fs::is_regular_file(filename));
    } catch (std::exception& e) {
        BOOST_THROW_EXCEPTION(e);
    }
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(createNoisedPatchtoSphereXml_test, DmOutputTestFixture)
{
    std::cout << "-- DmOutput: createNoisedPatchtoSphereXml_test" << std::endl;
    try
    {
        fs::path fitsFile{ "testFile.fits" };
        const fs::path filename = test_path / "NoisedPatchtoSphere.xml";
        std::string product_type = "DpdTwoDMassConvergencePatchesToSphere";
        int NResamples = 0;

        auto product = initProduct<dpdTwoDMassConvergencePatchesToSphere,
                                   twoDMassCollectConvergencePatchesToSphere,
                                   int>(product_type, NResamples);

        dm.createNoisedPatchtoSphereXml(product, fitsFile);
        writeProduct<dpdTwoDMassConvergencePatchesToSphere>(product, filename);
        BOOST_CHECK(fs::is_regular_file(filename));
    } catch (std::exception& e) {
        BOOST_THROW_EXCEPTION(e);
    }
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(createDenoisedPatchtoSphereXml_test, DmOutputTestFixture)
{
    std::cout << "-- DmOutput: createDenoisedPatchtoSphereXml_test"
               << std::endl;
    try
    {
        fs::path fitsFile { "testFile.fits" };
        const fs::path filename = test_path / "DeNoisedPatchtoSphere.xml";
        std::string product_type = "DpdTwoDMassConvergencePatchesToSphere";
        int NResamples = 0;

        auto product = initProduct<dpdTwoDMassConvergencePatchesToSphere,
                                   twoDMassCollectConvergencePatchesToSphere,
                                   int>(product_type, NResamples);

        dm.createDenoisedPatchtoSphereXml(product, fitsFile);
        writeProduct<dpdTwoDMassConvergencePatchesToSphere>(product, filename);
        BOOST_CHECK(fs::is_regular_file(filename));
    } catch (std::exception& e) {
        BOOST_THROW_EXCEPTION(e);
    }
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(createPatchtoSphereGalCountXml_test, DmOutputTestFixture)
{
    std::cout << "-- DmOutput: createPatchtoSphereGalCountXml_test"
              << std::endl;
    try
    {
        const fs::path filename = test_path / "GalCountPatchtoSphere.xml";
        fs::path fitsFile{ "testFile.fits" };
        std::string product_type = "DpdTwoDMassConvergencePatchesToSphere";
        int NResamples = 0;

        auto product = initProduct<dpdTwoDMassConvergencePatchesToSphere,
                                   twoDMassCollectConvergencePatchesToSphere,
                                   int>(product_type, NResamples);

        dm.createPatchtoSphereGalCountXml(product, fitsFile);
        writeProduct<dpdTwoDMassConvergencePatchesToSphere>(product, filename);
        BOOST_CHECK(fs::is_regular_file(filename));
    } catch (std::exception& e) {
        BOOST_THROW_EXCEPTION(e);
    }
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(createSNRPatchtoSphereOutputXml_test, DmOutputTestFixture)
{
    std::cout << "-- DmOutput: createSNRPatchtoSphereOutputXml_test"
              << std::endl;
    try
    {
        fs::path fitsFile{ "testFile.fits" };
        const fs::path filename = test_path / "SNRPatchtoSphere.xml";
        std::string product_type = "DpdTwoDMassConvergencePatchesToSphere";
        int NResamples = 2;
        auto product = initProduct<dpdTwoDMassConvergencePatchesToSphere,
                                   twoDMassCollectConvergencePatchesToSphere,
                                   int>(product_type, NResamples);

        dm.createSNRPatchtoSphereOutputXml(product, fitsFile);
        writeProduct<dpdTwoDMassConvergencePatchesToSphere>(product, filename);
        BOOST_CHECK(fs::is_regular_file(filename));
    } catch (std::exception& e) {
        BOOST_THROW_EXCEPTION(e);
    }
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(createPatchtoSphereProjCenterPosXml_test, DmOutputTestFixture)
{
    std::cout << "-- DmOutput: createPatchtoSphereProjCenterPosXml_test"
              << std::endl;
    try
    {
        fs::path fitsFile{ "testFile.fits" };
        const fs::path filename = test_path / "ProjCenterPosPatchtoSphere.xml";
        std::string product_type = "DpdTwoDMassConvergencePatchesToSphere";
        int NResamples = 0;

        auto product = initProduct<dpdTwoDMassConvergencePatchesToSphere,
                                   twoDMassCollectConvergencePatchesToSphere,
                                   int>(product_type, NResamples);

        dm.createPatchtoSphereProjCPosXml(product, fitsFile);
        writeProduct<dpdTwoDMassConvergencePatchesToSphere>(product, filename);
        BOOST_CHECK(fs::is_regular_file(filename));
    } catch (std::exception& e) {
        BOOST_THROW_EXCEPTION(e);
    }
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(createPatchtoSphereXml_test, DmOutputTestFixture)
{
    std::cout << "-- DmOutput: createPatchtoSphereXml_test" << std::endl;
    try
    {
        fs::path fitsFile1{ "testFile1.fits" };
        fs::path fitsFile2{ "testFile2.fits" };
        fs::path fitsFile3{ "testFile3.fits" };
        fs::path fitsFile4{ "testFile4.fits" };
        fs::path fitsFile5{ "testFile5.fits" };
        fs::path fitsFile6{ "testFile6.fits" };
        fs::path fitsFile7{ "testFile7.fits" };
        const fs::path filename = test_path / "PatchtoSphere.xml";
        std::string product_type = "DpdTwoDMassConvergencePatchesToSphere";
        int NResamples = 2;
        auto product = initProduct<dpdTwoDMassConvergencePatchesToSphere,
                                   twoDMassCollectConvergencePatchesToSphere,
                                   int>(product_type, NResamples);

        dm.createNoisedPatchtoSphereXml(product, fitsFile1);
        dm.createDenoisedPatchtoSphereXml(product, fitsFile2);
        dm.createPatchtoSphereGalCountXml(product, fitsFile3);
        dm.createSNRPatchtoSphereOutputXml(product, fitsFile4);
        dm.createPatchtoSphereProjCPosXml(product, fitsFile5);
        dm.createPatchtoSphereMCXml(product, fitsFile6);
        dm.createPatchtoSphereMCXml(product, fitsFile7);
        writeProduct<dpdTwoDMassConvergencePatchesToSphere>(product, filename);
        BOOST_CHECK(fs::is_regular_file(filename));
    } catch (std::exception& e) {
        BOOST_THROW_EXCEPTION(e);
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()
