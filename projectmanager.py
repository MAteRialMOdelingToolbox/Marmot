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
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--status', nargs='?',  default=False)
    # parser.add_argument('--branch', nargs='?',  default=False)
    # parser.add_argument('--checkout', nargs='?',  default=False)
    # parser.add_argument('--add', nargs='?',  default=False)
    # parser.add_argument('--push', nargs='?',  default=False)
    # parser.add_argument('--pull', nargs='?',  default=False)
    # parser.add_argument('--commit', nargs='?',  default=False)
    # args = parser.parse_args()



    directories = ['.']
    directories += list ( walklevel( './modules/materials/' , 1) )
    directories += list ( walklevel( './modules/elements/' , 1) )
    directories += ['./modules/bftMechanics'] 

    for d in directories:

        print(d)
        execs = 'bash -c "cd {:} && {:} " '.format(d, ' '.join ( sys.argv[1:] ) )
        os.system( execs )
            # if args.status:
                # os.system('bash -c "cd {:} && git status" '.format(d))

            # if args.branch:
                # os.system('bash -c "cd {:} && git branch {:}" '.format(d, args.branch))

            # if args.checkout:
                # os.system('bash -c "cd {:} && git checkout {:}" '.format(d, args.checkout))
            
            # if args.add:
                # os.system('bash -c "cd {:} && git add -u" '.format(d))

            # if args.commit:
                # os.system('bash -c "cd {:} && git commit -m \'{:}\'" '.format(d, args.commit))

            # if args.push:
                # os.system('bash -c "cd {:} && git push origin {:}" '.format(d, args.push))

            # if args.pull:
                # os.system('bash -c "cd {:} && git pull origin {:}" '.format(d, args.push))

