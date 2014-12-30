#! /usr/bin/env python

# This script is built to perform a pairwise comparison of samples
# containing locus information from DNA extraction
# The purpose is to score all samples against eachother to allow for
# easier diciphering of resampled individuals
# It then creates a series of graphics to help understand the data
from difflib import Differ
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from collections import Counter
from matplotlib.backends.backend_pdf import PdfPages
import math
import scipy as sp
import scipy.stats
import warnings
import argparse
import sys
import ast

warnings.simplefilter("error")
pd.options.mode.chained_assignment = None  # default='warn'

# Definitions
def parse_file(data):
	# map loci names to all alleles for each sample
	# return dictionary of data
	column_names = data[0].strip('\n').split(',')
	sample_dictionary = {}
	sample_id_list = []
	for index in range(1,len(data)):
		split_line = data[index].strip('\n').split(',')
		sample_id = split_line[0]
		sample_id_list.append(sample_id)
		locus_dictionary = {}
		for split_line_index in range(1,len(split_line)):
			allele_id = split_line[split_line_index]
			locus_id = column_names[split_line_index]
			if locus_id in locus_dictionary.keys():
				locus_dictionary[locus_id].append(allele_id)
			else:
				locus_dictionary[locus_id] = [allele_id]
		sample_dictionary[sample_id] = locus_dictionary
	return sample_dictionary

def make_pairwise(sample_dictionary):
	# create a list of all the possible comparisons of the samples
	# return list
	samples_pairwise = []
	for sample_index in range(0,len(sample_dictionary.keys())):
		sample_id = sample_dictionary.keys()[sample_index]
		for sample_index_2 in range(sample_index + 1,len(sample_dictionary.keys())):
			sample_id_2 = sample_dictionary.keys()[sample_index_2]
			samples_pairwise.append((sample_id,sample_id_2))
	return samples_pairwise

def score_pairs(samples_pairwise, sample_dictionary):
	# build scoring matrix from pairwise list
	# return list of scores
	score_data = []
	for samples in samples_pairwise:
		if (samples[0].find('virtsim') != -1 and str(samples[0].split('_')[-1]) == '0'
			) or (samples[0].find('virtsim') == -1
			):
			locus_dictionary_x = sample_dictionary[samples[0]]
			locus_dictionary_y = sample_dictionary[samples[1]]
			locus_pairwise = [list(Differ().compare(locus_dictionary_x[x],locus_dictionary_y[y]))
								for x in locus_dictionary_x.keys()
								for y in locus_dictionary_y.keys() if x == y]
			score_list = []
			for locus_list in locus_pairwise:
				use_list = bool([x for x in locus_list
							if x.split(' ')[-1] == '0'])
				if use_list != True:
					locus_scores = [x for x in locus_list
									 if ('+' in x or '-' in x)]
					score_list.append(len(locus_scores)/2)

			raw_score = sum(score_list)
			corrected_score = round(float(sum(score_list))/float(len(score_list)),2)
			score_data.append((samples[0],samples[1],raw_score,corrected_score))
	return score_data

def main(simdir, simlist):
	# set paths
	script_path = os.path.dirname(os.path.realpath(__file__))
	sim_path = simdir

	# make list of sim files
	list_dir = os.listdir(sim_path)
	if simlist.find('*') != -1:
		sim_list = [x for x in list_dir if x.find('.csv') != -1 and x.find('_score.csv') == -1]
	else:
		sim_list = ast.literal_eval(simlist)
	#sim_list = ['pairwise_error.csv', 'pairwise_full_noLW20.csv']

	# input model data to be used for real data analysis if it exists
	model_file = 'model_stats.tsv'
	if model_file in list_dir:
		model_input = pd.read_csv(sim_path + model_file, sep='\t', header=0, index_col=None, lineterminator='\n')
	else:
		columns = ['file','y_dist_max',
			'x_dist_max','x_dist_min','x_dist_mean','lower_bound','y_res_max',
			'x_res_max','x_res_min','x_res_mean','upper_bound']
		model_input = pd.DataFrame(index=[0], columns=columns)
		model_input = model_input.fillna(0)
	sim_model_data = []

	# start main processing loop for each sim file
	for file in sim_list:
		print file
		file = os.path.join(sim_path, file)
		# open sim file
		with open(file,'r') as f:
			data = f.readlines()

		# map loci names to all alleles for each sample
		sample_dictionary = parse_file(data)
		# create a list of all the possible comparisons of the samples
		samples_pairwise = make_pairwise(sample_dictionary)
		# build scoring matrix from pairwise list
		score_data = score_pairs(samples_pairwise, sample_dictionary)
		
		# create dataframe of scorring matrix
		score_df = pd.DataFrame(score_data,columns=['sample_0','sample_1','raw_score','corrected_score'])
		# create subset dataframe for use in creating histograms
		#hist_series = pd.DataFrame(score_df.corrected_score).convert_objects(convert_numeric=True)
		hist_series = score_df.corrected_score
		# Normalize the histogram data
		hs_sub_mean =  [round(x - hist_series.mean(),2) for x in hist_series]
		norm_hist_series = (hs_sub_mean/(hist_series.max() 
			- hist_series.min())
			)
		score_df['normalized_score'] = norm_hist_series

		# check if the file being processed is a sim or not
		resampled = []
		sim_OR_not = False
		for index, row in score_df.iterrows():
			individual_0 = '_'.join(row['sample_0'].split('_')[:-1])
			individual_1 = '_'.join(row['sample_1'].split('_')[:-1])
			# check to see if the file is a sim file or real data
			if individual_0.find('virtsim') != -1:
				sim_OR_not = True
			# if it is real data use the model file to check if there is a possible resample
			# a model file must exist for this to work properly
			if sim_OR_not == False:
				resampled.append(True if score_df.xs(index)['normalized_score'] < 
					model_input.lower_bound.mean() else False
					)
			# else it is a sims file and is checked for resamples
			else:
				resampled.append(True if individual_0 == individual_1 else False)
		# add the boolean resample list to the scoring dataframe
		score_df['resampled'] = resampled
		score_df = score_df.fillna('##')
		# output the scoring dataframe for safe keeping
		score_df.to_csv(file.rsplit('.', 1)[0] + '_score.csv')
		# make histogram of distinct individuals
		sim_stats = []
		distinct_individual_df = score_df[(score_df.resampled == False)]
		group_dist_in = pd.DataFrame(distinct_individual_df.groupby(['normalized_score']
			)['normalized_score'].count())
		group_dist_in.columns = ['count']
		group_dist_in = group_dist_in.reset_index()
		fig, ax = plt.subplots()
		# get stats for dist data
		y_dist_max = group_dist_in['count'].max()
		x_dist_max = group_dist_in['normalized_score'].max()
		x_dist_min = group_dist_in['normalized_score'].min()
		x_dist_mean = np.mean(group_dist_in['normalized_score'])
		x_dist_se = scipy.stats.sem(group_dist_in['normalized_score'])
		x_dist_count = group_dist_in['normalized_score'].count()
		n = len(group_dist_in['normalized_score'])
		confidence = 0.95
		x_dist_95_conf = x_dist_se * sp.stats.t._ppf((1 + confidence)/2., n - 1)
		lower_bound = x_dist_min - x_dist_95_conf
		# work with the simulated data
		if sim_OR_not == True:
			sim_stats = [file.split('/')[-1],str(y_dist_max),
				str(x_dist_max),str(x_dist_min),str(x_dist_mean),
				str(lower_bound)
				]
			bins = group_dist_in.normalized_score
			width = 0.7 * (bins[1] - bins[0])
			center = (bins[:-1] + bins[1:]) / 2
			plt.bar(group_dist_in.normalized_score, group_dist_in['count'], 
				align='center', width=width, color='b'
				)
			plt.plot([lower_bound,lower_bound], [y_dist_max,0], linestyle='dashed', color='b')

		# work with the Andrew's real data
		elif sim_OR_not == False:
			hist, bins = np.histogram(distinct_individual_df.normalized_score, bins=20)
			width = 0.7 * (bins[1] - bins[0])
			center = (bins[:-1] + bins[1:]) / 2
			colors = []
			for bar in center:
				if bar < model_input.lower_bound.mean():
					color = 'r'
				else:
					color = 'b'
				colors.append(color)
			plt.bar(center, hist, align='center',
				width=width, color=colors
				)
			bins = group_dist_in['normalized_score']
			# determine color of shaded area
			if model_input.upper_bound.mean() <= model_input.lower_bound.mean():
				facecolor = 'g'
			if model_input.upper_bound.mean() > model_input.lower_bound.mean():
				facecolor = 'y'
			# plot a shaded area for the average max resampled and average min distinct from model
			plt.axvspan(model_input.upper_bound.mean(),
				model_input.lower_bound.mean(), facecolor=facecolor, alpha=0.5)
			plt.plot([model_input.lower_bound.mean(),model_input.lower_bound.mean()],
				[y_dist_max,0], linestyle='dashed', color='b'
				)
			plt.plot([model_input.upper_bound.mean(),model_input.upper_bound.mean()],
			[y_dist_max,0], linestyle='dashed', color='r'
				)

		# make histogram of resampled individuals
		resampled_individual_df = score_df[(score_df.resampled == True)]
		# add general virtsim id for creation of second graph
		general_sample = resampled_individual_df['sample_0'].apply(lambda x: '_'.join(x.split('_')[0:2]))
		resampled_individual_df['general_sample'] = general_sample
		found_resamples = resampled_individual_df.shape[0]
		if found_resamples != 0:
			group_res_in = pd.DataFrame(resampled_individual_df.groupby(['normalized_score']
				)['normalized_score'].count())
			group_res_in.columns = ['count']
			group_res_in = group_res_in.reset_index()

			quant_res_in = pd.DataFrame(resampled_individual_df.groupby(['general_sample']
				)['general_sample'].count())
			quant_res_in.columns = ['count']
			quant_res_in = quant_res_in.reset_index()
			# get stats for res data
			y_res_max = group_res_in['count'].max()
			x_res_max = group_res_in['normalized_score'].max()
			x_res_min = group_res_in['normalized_score'].min()
			x_res_mean = group_res_in['normalized_score'].mean()
			if sim_OR_not == True:
				x_res_se = scipy.stats.sem(group_res_in['normalized_score'])
				if math.isnan(x_res_se) == True:
					x_res_std = 0
				x_res_count = group_res_in['normalized_score'].count()
				n = len(group_res_in['normalized_score'])
				confidence = 0.95
				x_res_95_conf = x_res_se * sp.stats.t._ppf((1 + confidence)/2., n - 1)
				upper_bound = x_res_max + x_res_95_conf
				sim_stats = sim_stats + [str(y_res_max),str(x_res_max),str(x_res_min),
					str(x_res_mean),str(upper_bound)
					]
			bins = group_res_in['normalized_score']
			# plot the histogram data for resampled
			plt.bar(group_res_in.normalized_score, group_res_in['count'],
				align='center', width=width, color='r'
				)
			# determine color of shaded area
			if sim_OR_not == True:
				# plot dashed line for max resampled bin
				plt.plot([upper_bound,upper_bound], [y_dist_max,0],
					linestyle='dashed', color='r'
					)
				if upper_bound <= lower_bound:
					facecolor = 'g'
				if upper_bound > lower_bound:
					facecolor = 'y'
				# plot a shaded area for the max resampled and min distinct
				plt.axvspan(upper_bound, lower_bound, facecolor=facecolor, alpha=0.5)
			
		sim_model_data.append(sim_stats)

		plt.xlim([-1,1])
		plt.xlabel('score bins')
		plt.ylabel('frequency')
		plt.title(file.split('.')[0].split('/')[-1])
		pp = PdfPages(file.rsplit('.',1)[0] + '.pdf')
		pp.savefig(fig)
		plt.clf()

		if found_resamples != 0:
			fig, ax = plt.subplots()
			quant_res_in = pd.DataFrame(resampled_individual_df.groupby(['general_sample']
				)['general_sample'].count())
			quant_res_in.columns = ['count']
			quant_res_in = quant_res_in.reset_index()
			resample_list = [str(x[0]) + ' ' + str(x[1]) 
				for x in zip(quant_res_in.general_sample,quant_res_in['count'])
				]
			quant_res_in.columns = ['num_indv','num_resamp']
			quant_double_group = quant_res_in.groupby(['num_resamp']
				)['num_indv'].count().reset_index()
			quant_double_group.columns = ['num_resamp', 'num_indv']
			x_max = quant_double_group['num_resamp'].max()
			x_min = quant_double_group['num_resamp'].min()
			y_max = quant_double_group['num_indv'].max()
			bins = quant_double_group['num_resamp']
			# plot the histogram data for resampled
			plt.bar(bins, quant_double_group['num_indv'],
				align='center', width=width, color='r'
				)
			plt.xticks(np.arange(0,x_max + 2,1))
			plt.yticks(np.arange(0,y_max + 2,1))
			plt.xlabel('# of resamples/individual')
			plt.ylabel('# of individuals resampled')
			ax.text(x_max + 0.25, 0.25, '\n'.join(resample_list), style='italic',
			    bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
			pp.savefig(fig)
			plt.clf()
		pp.close()
		plt.close('all')

	write2model = False
	if write2model == True:
		with open(sim_path + 'model_stats.tsv','w') as o:
			o.write('\t'.join(['file','y_dist_max',
				'x_dist_max','x_dist_min','x_dist_mean','lower_bound','y_res_max',
				'x_res_max','x_res_min','x_res_mean','upper_bound']
				) + '\n')
			for stats in sim_model_data:
				o.write('\t'.join(stats) + '\n')


# Start Main
if __name__ == '__main__':
# collect arguments from commandline
	parser = argparse.ArgumentParser(description='insert program description')
	parser.add_argument('-s','--simdir', help='path to simulation file directory.', 
		required=False, default=os.path.dirname(os.path.realpath(__file__))
		)
	parser.add_argument('-l','--simlist', help='list of simulation files.', 
		required=True
		)
	args = vars(parser.parse_args())
	# check to make sure that enough arguments were passed before proceeding
	if len(sys.argv) < 2:
		sys.exit("Missing flags, type \"--help\" for assistance...")
	main(args['simdir'], args['simlist'])
