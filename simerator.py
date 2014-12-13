import sys
from subprocess import Popen, PIPE
import argparse


def main(simtype, infile, boutlimit,
	samplelimit, popsize, outfile, iterations):
	for i in range(0,int(iterations)):
		print i
		run = Popen(['python', 'BeanBag.py', '--simtype',simtype,'--infile',
			infile,'--popsize',popsize,'--samplelimit',samplelimit,'--boutlimit',boutlimit,'--outfile',
			outfile + '_I' + str(i) + '.csv'], stdout=PIPE)
		print run.stdout.read()


#Start Main
if __name__ == '__main__':
	# collect arguments from commandline
	parser = argparse.ArgumentParser(description='insert program description')
	parser.add_argument('--simtype', help='Specify which type of simulation type; virtpop, simfile, genepop.',
		required=True, default=False
		)
	parser.add_argument('--infile',	help='Specify path to input file; GENEPOP, SIMFILE, or no input.',
		required=True, default=False
		)
	parser.add_argument('--boutlimit',help='Total number of sampling bouts for a season.',
		required=True, default=False
		)
	parser.add_argument('--samplelimit', help='Total number of samples for a season.',
		required=True, default=False
		)
	parser.add_argument('--popsize', help='Specify the size of a virtual population.',
		required=True, default=False
		)
	parser.add_argument('--outfile', help='Specify the output file.',
		required=True, default=False
		)
	parser.add_argument('--iterations',	help='how many iterations.',
		required=True, default=False
		)
	args = vars(parser.parse_args())
	# check to make sure that enough arguments were passed before proceeding
	if len(sys.argv) < 2:
		sys.exit("Missing flags, type \"--help\" for assistance...")
	main(args['simtype'], args['infile'], args['boutlimit'], args['samplelimit'],
		args['popsize'], args['outfile'], args['iterations']
		)
