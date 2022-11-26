/**
 * @file src/program/CartesianWL.cpp
 * @date 02/14/22
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

#include "ElementsKernel/ProgramHeaders.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianProcessor.h"
#include <chrono>
#include <map>
#include <string>
#include <utility>

using LE3_2D_MASS_WL_CARTESIAN::CartesianProcessor;

class CartesianWL: public Elements::Program
{

public:

    std::pair<OptionsDescription, PositionalOptionsDescription> defineProgramArguments()
            override
    {
        OptionsDescription desc { };
        auto add = desc.add_options();
        add("paramfile", value<std::string>()->required(), "Input parameters file");
        add("workdir", value<std::string>()->required(), "Working directory");
        add("shear", value<std::string>(), "Input shear catalog");
        add("clusterCatalog", value<std::string>(), "Input cluster catalog");
        add("mc", "Do produce Monte Carlo maps");

        PositionalOptionsDescription pos_desc;
        pos_desc.add("paramfile", 1);
        return std::make_pair(desc, pos_desc);
    }

    ExitCode mainMethod(std::map<std::string, VariableValue>& args) override
    {
        ////////////////////////////////////////////////////////////////////////
        // Start of program
        ////////////////////////////////////////////////////////////////////////
        Logging logger = Logging::getLogger("CartesianWL");
        auto start_time = std::chrono::system_clock::now();
        logger.info("#");
        logger.info("# Entering mainMethod()");
        logger.info("#");

        ////////////////////////////////////////////////////////////////////////
        // Run process
        ////////////////////////////////////////////////////////////////////////

        try
        {
            CartesianProcessor cp(args);
            cp.process();
        }
        catch(ExitCode& exitCode)
        {
            int code = static_cast<int>(exitCode);
            logger.fatal("Exit code: " + std::to_string(code));
            return exitCode;
        }

        ////////////////////////////////////////////////////////////////////////
        // End of program
        ////////////////////////////////////////////////////////////////////////
        auto end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_time - start_time;
        logger.info() << "Elapsed time: " << elapsed_seconds.count() << "s";
        logger.info("#");
        logger.info("# Exiting mainMethod()");
        logger.info("#");

        return ExitCode::OK;
    }

};

MAIN_FOR(CartesianWL)

