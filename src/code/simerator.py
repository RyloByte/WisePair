import sys
from subprocess import Popen, PIPE
import argparse


def main(simtype, infile, errfile, boutlimit,
         samplelimit, popsize, perpop, perbout, outfile, iterations):
    for i in range(0, int(iterations)):
        print i
        run = Popen(['python', 'BeanBag.py', '-t', simtype, '-i', infile, '-e', errfile,
                     '-p', popsize, '-r', perpop, '-l', str(perbout), '-s', samplelimit, '-b', boutlimit, '-o',
                     outfile + '_I' + str(i) + '.csv'], stdout=PIPE)
        print run.stdout.read()


# Start Main
if __name__ == '__main__':
    # collect arguments from commandline
    parser = argparse.ArgumentParser(description='insert program description')
    parser.add_argument('-t', '--simtype', help='Specify which type of simulation type; virtpop, simfile, genepop.',
                        required=True
                        )
    parser.add_argument('-i', '--infile', help='Specify path to input file; GENEPOP, SIMFILE, or no input.',
                        required=True
                        )
    parser.add_argument('-e', '--errfile', help='Specify path to input error file; JSON format',
                        required=True
                        )
    parser.add_argument('-b', '--boutlimit', help='Total number of sampling bouts for a season.',
                        required=True
                        )
    parser.add_argument('-s', '--samplelimit', help='Total number of samples for a season.',
                        required=True
                        )
    parser.add_argument('-p', '--popsize', help='Specify the size of a virtual population.',
                        required=True
                        )
    parser.add_argument('-r', '--perpop', help='Specify the percent of a population present at bout.',
                        required=False, default=1
                        )
    parser.add_argument('-l', '--perbout', help='Number of samples for each bout.',
                        required=False, default=False
                        )
    parser.add_argument('-o', '--outfile', help='Specify the output file.',
                        required=True
                        )
    parser.add_argument('-m', '--iterations', help='Number of iterations to perform.',
                        required=True
                        )
    args = vars(parser.parse_args())
    # check to make sure that enough arguments were passed before proceeding
    if len(sys.argv) < 7:
        sys.exit("Missing flags, type \"--help\" for assistance...")
    main(args['simtype'], args['infile'], args['errfile'], args['boutlimit'], args['samplelimit'],
         args['popsize'], args['perpop'], args['perbout'], args['outfile'], args['iterations']
         )
