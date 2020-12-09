import os
import sys
import subprocess

def walk_modules(rootdir, levels=1):

    rootdir = rootdir.rstrip(os.path.sep)
    assert os.path.isdir(rootdir)
    num_sep = rootdir.count(os.path.sep)

    modules = []

    for directory, subdirs, files in os.walk(rootdir):

        num_sep_this = directory.count(os.path.sep)
        if num_sep + levels <= num_sep_this:
            continue

        if 'module.cmake' in files:
            modules.append (directory)
                
    return modules

if __name__ == '__main__':

    projects = ['.'] 
    projects += walk_modules( './modules/', 3 )
    for proj in projects:
        print(proj)
        p = subprocess.Popen(sys.argv[1:], cwd=proj)
        p.wait()

