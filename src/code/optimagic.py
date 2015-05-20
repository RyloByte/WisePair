#! /usr/bin/env python

#############################
# Import specific libraries #
#############################

import sys
import argparse
from os.path import dirname, realpath, join, exists
from os import listdir, makedirs
import pandas as pd
from subprocess import Popen, PIPE
from itertools import product


def count_resampled(simdir_list, resampmin, nummin):
    stats_df = pd.DataFrame(columns=['tot_resamp','mean_rs_per_ind','met_resamp_min','mean_rs_per_mrm'])
    for file in simdir_list:
        score_df = pd.read_csv(file)
        resampled_df = score_df[['sample_0','sample_1','normalized_score','resampled']
                ].query('resampled == True'
                )
        sample_0_list = [x.rsplit('_', 1)[0] for x in resampled_df.sample_0]
        sample_1_list = [x.rsplit('_', 1)[0] for x in resampled_df.sample_1]
        count_tup_list = [(x, sample_1_list.count(x)) for x in set(sample_0_list)
            ]
        if count_tup_list != []:
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

def run_sim(simdir, simtype, simfile, errfile, popsize, perpop, samplim, boutlim, iterlim):
    print 'Using Simerator to run %s iterations...' % str(iterlim)
    print 'Running BeanBag simulations...'
    om_work_dir = set_simdir(simdir, 'OM_simdir')
    outfile = join(om_work_dir, ''.join(['P',popsize,'_S',samplim,'_B',boutlim,'_R', perpop]))
    run = Popen(['python', 'simerator.py', '-t',simtype,'-i',simfile,
            '-e', errfile,'-p',popsize, '-r',perpop,'-s',samplim,'-b',boutlim,
            '-m',iterlim,'-o',outfile], stdout=PIPE)
    run.stdout.read()
    return om_work_dir

def run_wise(simdir, simlist):
    print 'Running WisePair scoring algorithm...\n'
    run = Popen(['python', 'wise_pair.py', '-s', simdir, '-l', str(simlist)], stdout=PIPE)
    run.stdout.read()
    
def optimaker(simdir, popsize, perpop, samplim, boutlim, resampmin, nummin, 
    runsim, iterlim, simtype, simfile, errfile
    ):
    if runsim == 'True':
        om_work_dir = run_sim(simdir, simtype, simfile, errfile, 
            popsize, perpop, samplim, boutlim, iterlim
            )
        simdir = om_work_dir
        simdir_list = [join(simdir, x) for x in listdir(simdir) if x.find('_score.csv') == -1 and
            x.find('.pdf') == -1 and    x.find('P' + str(popsize)) != -1 and
            x.find('S' + str(samplim)) != -1 and x.find('B' + str(boutlim)) != -1
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

    print '#########################################################################################'
    print 'Criteria:\n\t%s Resampled Individutals\n\t%s Times Resamples' % (resampmin, nummin)
    print '\t%s Number of Bouts\n\t%s Samples per Bout \n\tPopulation of %s\n' % (boutlim, str(int(samplim)/int(boutlim)), popsize)
    #print find_resampled[0]
    print find_resampled[1]
    print '#########################################################################################'
    print '\n\n'
    return find_resampled[1]

def main(simdir, popsize, perpop, sampran, boutran, resampmin,
    numminran, runsim, iterlim, simtype, simfile, errfile
    ):
    if sampran.find(',') != -1:
        sampran_list = range(int(sampran.split(',')[0]), int(sampran.split(',')[1]) + 1)
    else:
        sampran_list = [sampran]
    if boutran.find(',') != -1:
        boutran_list = range(int(boutran.split(',')[0]), int(boutran.split(',')[1]) + 1)
    else:
        boutran_list = [boutran]
    if numminran.find(',') != -1:
        numminran_list = range(int(numminran.split(',')[0]), int(numminran.split(',')[1]) + 1)
    else:
        numminran_list = [numminran]

    range_tuples = list(product(sampran_list,boutran_list,numminran_list))
    print 'Running %s simulations.' % (str(len(range_tuples)))
    good_sim_stats_df = pd.DataFrame(columns=['population', 'number_of_samples',
        'number_of_bout', 'min_resampled', 'min_times_resampled',
        'tot_resamp', 'mean_rs_per_ind', 'met_resamp_min', 'mean_rs_per_mrm']
        )
    other_sim_stats_df = pd.DataFrame(columns=['population', 'number_of_samples',
        'number_of_bout', 'min_resampled', 'min_times_resampled',
        'tot_resamp', 'mean_rs_per_ind', 'met_resamp_min', 'mean_rs_per_mrm']
        )
    for range_tuple in range_tuples:
        boutlim = str(range_tuple[1])
        samplim = str(int(range_tuple[0]) * int(boutlim))
        nummin = str(range_tuple[2])
        sim_stats_mean = optimaker(simdir, popsize, perpop, samplim, boutlim, resampmin, nummin, 
            runsim, iterlim, simtype, simfile, errfile)
        if (sim_stats_mean[2] >= int(resampmin) and sim_stats_mean[3] >= int(nummin)):
            good_sim_stats_df.loc[len(good_sim_stats_df) + 1] = [popsize, samplim, boutlim, resampmin, nummin,
                round(sim_stats_mean[0],2), round(sim_stats_mean[1],2), round(sim_stats_mean[2],2),
                round(sim_stats_mean[3],2)
                ]
        else:
            other_sim_stats_df.loc[len(other_sim_stats_df) + 1] = [popsize, samplim, boutlim, resampmin, nummin,
                round(sim_stats_mean[0],2), round(sim_stats_mean[1],2), round(sim_stats_mean[2],2),
                round(sim_stats_mean[3],2)
                ]
    good_sim_stats_df.to_csv(simdir + '/optimagic_output.csv')
    other_sim_stats_df.to_csv(simdir + '/optimagic_output_BAD.csv')



if __name__ == '__main__':
    # collect arguments from commandline
    parser = argparse.ArgumentParser(description='insert program description')
    parser.add_argument('-s','--simdir', help='path to simulation file directory.', 
        required=False, default=dirname(realpath(__file__))
        )
    parser.add_argument('-p','--popsize', help='population size.', 
        required=True
        )
    parser.add_argument('-x','--perpop', help='population size present at site.', 
        required=False, default=100
        )
    parser.add_argument('-l','--sampran', help='range of samples per bout.', 
        required=True
        )
    parser.add_argument('-b','--boutran', help='range of bouts.', 
        required=True
        )
    parser.add_argument('-r','--resampmin', help='minimum individuals resampled over season.', 
        required=True
        )
    parser.add_argument('-n','--numminran', help='range for minimum number of resamples per resampled individual.', 
        required=True
        )
    parser.add_argument('-u','--runsim', help='use simerator to run simulations.', 
        required=False, default=False
        )
    parser.add_argument('-m','--iterlim', help='number of iterations during simulation mode.', 
        required=False, default=False
        )
    parser.add_argument('-t','--simtype', help='specitfy type of input simulation format: virtpop, simfile, genepop.', 
        required=False, default=False
        )
    parser.add_argument('-i','--simfile', help='path to simulation input file.', 
        required=False, default=False
        )
    parser.add_argument('-e','--errfile', help='Specify path to input error file; JSON format', 
        required=False, default=False
        )
    args = vars(parser.parse_args())
    # check to make sure that enough arguments were passed before proceeding
    if len(sys.argv) < 2:
        sys.exit("Missing flags, type \"--help\" for assistance...")
    main(args['simdir'], args['popsize'], args['perpop'], args['sampran'], args['boutran'],
        args['resampmin'], args['numminran'], args['runsim'], args['iterlim'],
        args['simtype'], args['simfile'], args['errfile']
        )
