/**
 * @file LE3_2D_MASS_WL_CARTESIAN/CartesianProcessor.h
 * @date 02/15/22
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

#ifndef _LE3_2D_MASS_WL_CARTESIAN_CARTESIANPROCESSOR_H
#define _LE3_2D_MASS_WL_CARTESIAN_CARTESIANPROCESSOR_H

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include "ElementsKernel/ProgramHeaders.h"
#include <string>
#include <map>

using Elements::ExitCode;
using boost::filesystem::path;
using boost::program_options::value;
using boost::program_options::variable_value;

namespace LE3_2D_MASS_WL_CARTESIAN
{

/**
 * @class CartesianProcessor
 * @brief
 *
 */

class CartesianProcessor
{

public:

    /**
     * @brief Constructor
     */
    explicit CartesianProcessor(std::map<std::string, variable_value>& args);

    /**
     * @brief Constructor
     */
    explicit CartesianProcessor(std::map<std::string, std::string>& args);

    /**
     * @brief Parse the options for the parameter type
     */
    void parseOptions();

    /**
     * @brief Process the operation
     */
    void process();

private:
    std::map<std::string, std::string> m_args;
    path m_workdir;
    path m_paramfile;
    std::string m_paramtype;
    path m_shearCatalogPath;
    path m_shearMapPath;
    path m_clusterCatalog;
    bool m_produceMcMaps;


};
// End of CartesianProcessor class

}// namespace LE3_2D_MASS_WL_CARTESIAN

#endif // _LE3_2D_MASS_WL_CARTESIAN_CARTESIANPROCESSOR_H
