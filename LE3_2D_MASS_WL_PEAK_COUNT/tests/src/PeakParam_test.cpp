/*
 * Copyright (C) 2012-2020 Euclid Science Ground Segment
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
 */

/**
 * @file tests/src/PeakParam_test.cpp
 * @date 09/17/19
 * @author user
 */

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>

#include "LE3_2D_MASS_WL_PEAK_COUNT/PeakParam.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

#include "ElementsServices/DataSync.h"

#include <ios>
#include <sstream>
#include <iostream>

using namespace LE3_2D_MASS_WL_PEAK_COUNT;

using namespace ElementsServices::DataSync;

using boost::filesystem::path;

//-----------------------------------------------------------------------------

struct PeakParamDataSyncFixture
{
    DataSync sync;
    path xmlPConvFileName, xmlMAFileName;
    PeakParamDataSyncFixture() :
            sync("LE3_2D_MASS_WL_PEAK_COUNT/datasync_webdav.conf",
                 "LE3_2D_MASS_WL_PEAK_COUNT/test_file_list.txt"),
            xmlPConvFileName(sync.absolutePath("PConvparam_test.xml")),
            xmlMAFileName(sync.absolutePath("PMassAp_test.xml"))
    {
        sync.download();
    }
};

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (PeakParam_test, PeakParamDataSyncFixture)
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(xmlPConvFile_test, PeakParamDataSyncFixture)
{
    std::cout << "-- PARAM: .xml file for Wavelet" << std::endl;
    PeakParam Param;
    Param.readWaveletPeakXML(xmlPConvFileName);
    BOOST_CHECK_EQUAL(Param.getNPeakScale(), 5);
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(xmlMAFile_test, PeakParamDataSyncFixture)
{
    std::cout << "-- PARAM: .xml file for Mass Aperture" << std::endl;
    PeakParam Param;
    Param.readMassApertureXML(xmlMAFileName);
    std::vector<double> apPeakRadius = Param.getApPeakRadius();
    BOOST_CHECK_EQUAL(apPeakRadius[0], 1);
    BOOST_CHECK_EQUAL(apPeakRadius[1], 2);
    BOOST_CHECK_EQUAL(apPeakRadius[2], 3);
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(readPeakParamFile_test, PeakParamDataSyncFixture)
{
    std::cout << "-- PARAM: readPeakParamFile_test" << std::endl;
    PeakParam Param;
    readPeakParamFile(xmlMAFileName, Param);
    std::vector<double> apPeakRadius = Param.getApPeakRadius();
    BOOST_CHECK_EQUAL(apPeakRadius[0], 1);
    BOOST_CHECK_EQUAL(apPeakRadius[1], 2);
    BOOST_CHECK_EQUAL(apPeakRadius[2], 3);
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(readPeakParamFile_test2, PeakParamDataSyncFixture)
{
    std::cout << "-- PARAM: readPeakParamFile_test2" << std::endl;
    PeakParam Param;
    readPeakParamFile(xmlPConvFileName, Param);
    BOOST_CHECK_EQUAL(Param.getNPeakScale(), 5);
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(readPeakParamFile_test3, PeakParamDataSyncFixture)
{
    std::cout << "-- PARAM: readPeakParamFile_test3" << std::endl;
    auto fakeFilePath = fs::path("wrong_path/wrong_file.ini");
    PeakParam Param;
    BOOST_CHECK_THROW(readPeakParamFile(fakeFilePath, Param),
                      Elements::Exception);
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(readWrongFile_test)
{
    PeakParam Param;
    Param.readWaveletPeakXML("wrong_file_path");
    Param.readMassApertureXML("wrong_file_path");
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()
