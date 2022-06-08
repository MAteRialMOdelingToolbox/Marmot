# /* ---------------------------------------------------------------------
# *                                       _
# *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
# * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
# * | | | | | | (_| | |  | | | | | | (_) | |_
# * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
# *
# * Unit of Strength of Materials and Structural Analysis
# * University of Innsbruck,
# * 2020 - today
# *
# * festigkeitslehre@uibk.ac.at
# *
# * Alexander Dummer alexander.dummer@uibk.ac.at
# *
# * This file is part of the MAteRialMOdellingToolbox (marmot).
# *
# * This library is free software; you can redistribute it and/or
# * modify it under the terms of the GNU Lesser General Public
# * License as published by the Free Software Foundation; either
# * version 2.1 of the License, or (at your option) any later version.
# *
# * The full text of the license can be found in the file LICENSE.md at
# * the top level directory of marmot.
# * ---------------------------------------------------------------------
# */

import os
import sys

from projectmanager import walk_modules

rootDir = os.getcwd()
preCommitInstallCommand = "cp {root}/.pre-commit-config.yaml .".format(root=rootDir)
preCommitInstallCommand += " && pre-commit install"
preCommitInstallCommand += " && rm .pre-commit-config.yaml"

if __name__ == "__main__":

    # check if pre-commit is installed already
    try:
        os.system("pre-commit --help >/dev/null 2>&1")
    except:
        os.system("pip install pre-commit")

    os.system("pre-commit install")
    projects = walk_modules("./modules/", 3)
    for proj in projects:
        print(proj)
        os.chdir(proj)
        os.system(preCommitInstallCommand)
        os.chdir(rootDir)
