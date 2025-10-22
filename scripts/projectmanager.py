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
# * Matthias Neuner matthias.neuner@uibk.ac.at
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
import subprocess
import sys


def walk_modules(rootdir, levels=1):

    rootdir = rootdir.rstrip(os.path.sep)
    assert os.path.isdir(rootdir)
    num_sep = rootdir.count(os.path.sep)

    modules = []

    for directory, subdirs, files in os.walk(rootdir):

        num_sep_this = directory.count(os.path.sep)
        if num_sep + levels <= num_sep_this:
            continue

        if "module.cmake" in files:
            modules.append(directory)

    return modules


if __name__ == "__main__":

    projects = ["."]
    projects += walk_modules("./modules/", 3)
    for proj in projects:
        print(proj)
        p = subprocess.Popen(sys.argv[1:], cwd=proj)
        p.wait()
