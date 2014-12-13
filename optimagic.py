#############################
# Import specific libraries #
#############################

import sys
import argparse
from os.path import dirname, realpath, join, exists
from os import listdir, makedirs
import pandas as pd
from subprocess import Popen, PIPE


def count_resampled(simdir_list, resampmin, nummin):
	stats_df = pd.DataFrame(columns=['tot_resamp','mean_rs_per_ind','resamp_min','mean_rs_per_rm'])
	for file in simdir_list:
		score_df = pd.read_csv(file)
		resampled_df = score_df[['sample_0','sample_1','normalized_score','resampled']
				].query('resampled == True'
				)
		sample_0_list = [x.rsplit('_', 1)[0] for x in resampled_df.sample_0]
		sample_1_list = [x.rsplit('_', 1)[0] for x in resampled_df.sample_1]
		count_tup_list = [(x, sample_1_list.count(x)) for x in set(sample_0_list)
			]
		resamp_count_df = pd.DataFrame(count_tup_list, 
			columns=['virtsim_id', 'resample_count']
			)
		test_resamples = resamp_count_df[resamp_count_df['resample_count'] >= int(nummin)]
		stats_df.loc[len(stats_df)] = [len(resamp_count_df), round(resamp_count_df.resample_count.mean(),1),
			len(test_resamples), round(test_resamples.resample_count.mean(),1)
			]
	stats_df.fillna(0, inplace=True)
	return stats_df, stats_df.mean()

def set_simdir(simdir, om_dir_name):
	om_work_dir = join(simdir, om_dir_name)
	if not exists(om_work_dir):
		makedirs(om_work_dir)
	return om_work_dir

def run_sim(simdir, simtype, simfile, popsize, samplim, boutlim, iterlim):
	print 'Running BeanBag simulations...\n'
	om_work_dir = set_simdir(simdir, 'OM_simdir')
	outfile = join(om_work_dir, ''.join(['P',popsize,'_S',samplim,'_B',boutlim]))
	run = Popen(['python', 'simerator.py', '--simtype',simtype,'--infile',
			simfile,'--popsize',popsize,'--samplelimit',samplim,'--boutlimit',boutlim,
			'--iterations',iterlim,'--outfile',outfile], stdout=PIPE)
	run.stdout.read()
	return om_work_dir

def run_wise(simdir, simlist):
	print 'Running WisePair scoring algorithm...\n'
	run = Popen(['python', 'wise_pair.py', '-s', simdir, '-l', str(simlist)], stdout=PIPE)
	run.stdout.read()
	
def main(simdir, popsize, samplim, boutlim, resampmin, nummin, 
	runsim, iterlim, simtype, simfile
	):
	if runsim == 'True':
		om_work_dir = run_sim(simdir, simtype, simfile, popsize, samplim, boutlim, iterlim)
		simdir = om_work_dir
		simdir_list = [join(simdir, x) for x in listdir(simdir) if x.find('_score.csv') == -1 and 
			x.find('P' + str(popsize)) != -1 and x.find('S' + str(samplim)) != -1 and 
			x.find('B' + str(boutlim)) != -1
			]
		run_wise(simdir, simdir_list)
	
	simdir_list = [join(simdir, x) for x in listdir(simdir) if x.find('_score.csv') != -1 and 
		x.find('P' + str(popsize)) != -1 and x.find('S' + str(samplim)) != -1 and 
		x.find('B' + str(boutlim)) != -1
		]

	if simdir_list == []:
		print 'No score files found in %s' % simdir
		sys.exit()
	find_resampled = count_resampled(simdir_list, resampmin, nummin)
	print find_resampled[0]
	print find_resampled[1]




if __name__ == '__main__':
	# collect arguments from commandline
	parser = argparse.ArgumentParser(description='insert program description')
	parser.add_argument('-s','--simdir', help='path to simulation file directory.', 
		required=False, default=dirname(realpath(__file__))
		)
	parser.add_argument('-p','--popsize', help='population size.', 
		required=True
		)
	parser.add_argument('-l','--samplim', help='total sample limit', 
		required=True
		)
	parser.add_argument('-b','--boutlim', help='bout limit.', 
		required=True
		)
	parser.add_argument('-r','--resampmin', help='minimum individuals resampled over season.', 
		required=True
		)
	parser.add_argument('-n','--nummin', help='minimum number of resamples per resampled individual.', 
		required=True
		)
	parser.add_argument('-u','--runsim', help='use simerator to run simulations.', 
		required=False, default=False
		)
	parser.add_argument('-i','--iterlim', help='number of iterations during simulation mode.', 
		required=False, default=False
		)
	parser.add_argument('-t','--simtype', help='specitfy type of input simulation format: virtpop, simfile, genepop.', 
		required=False, default=False
		)
	parser.add_argument('-f','--simfile', help='path to simulation input file.', 
		required=False, default=False
		)
	args = vars(parser.parse_args())
	# check to make sure that enough arguments were passed before proceeding
	if len(sys.argv) < 2:
		sys.exit("Missing flags, type \"--help\" for assistance...")
	main(args['simdir'], args['popsize'], args['samplim'], args['boutlim'],
		args['resampmin'], args['nummin'], args['runsim'], args['iterlim'],
		args['simtype'], args['simfile']
		)
