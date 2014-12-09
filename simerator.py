from sys import argv
from subprocess import Popen, PIPE
import click


@click.command()
@click.option('--simtype',
	help='Specify which type of simulation type; virtpop, simfile, genepop.'
	)
@click.option('--infile',
	help='Specify path to input file; GENEPOP, SIMFILE, or no input.'
	)
@click.option('--boutlimit',
	help='Total number of sampling bouts for a season.'
	)
@click.option('--samplelimit',
	help='Total number of samples for a season.'
	)
@click.option('--popsize',
	help='Specify the size of a virtual population.'
	)
@click.option('--outfile',
	help='Specify the output file.'
	)
def main(simtype, infile, boutlimit,
	samplelimit, popsize, outfile):
	for i in range(0,int(iterations)):
		print i
		run = subprocess.Popen(['python', 'BeanBag_click_it.py', '--simtype',simfile,'--infile',
			infile,'--popsize',popsize,'--samplelimit',samplelimit,'--boutlimit',boutlimit,'--outfile',
			outfile + str(i) + '.csv'], stdout=subprocess.PIPE)
		print run.stdout.read()