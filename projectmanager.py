import os
import argparse
import sys

def walklevel(some_dir, level=1):
    some_dir = some_dir.rstrip(os.path.sep)
    assert os.path.isdir(some_dir)
    num_sep = some_dir.count(os.path.sep)
    for root, dirs, files in os.walk(some_dir):
        # yield root, dirs, files
        dirs  = [ os.path.join(some_dir, d ) for d in dirs  if os.path.isdir(os.path.join(some_dir, d ) )  ]
        return dirs
        # yield r dirs, files
        num_sep_this = root.count(os.path.sep)
        if num_sep + level <= num_sep_this:
            del dirs[:]

if __name__ == '__main__':

    directories = ['.']
    directories += list ( walklevel( './modules/materials/' , 1) )
    directories += list ( walklevel( './modules/elements/' , 1) )
    directories += ['./modules/MarmotMechanics'] 
    directories += ['./modules/MarmotFiniteElementCore'] 
    directories += ['./modules/MarmotCosseratCore'] 
    directories += ['./modules/MarmotMicromorphicCore'] 

    for d in directories:

        print(d)
        execs = 'bash -c "cd {:} && {:} " '.format(d, ' '.join ( sys.argv[1:] ) )
        os.system( execs )
