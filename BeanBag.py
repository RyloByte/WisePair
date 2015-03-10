#! /usr/bin/env python

from random import randint, choice, sample
from collections import Counter
import sys
import math
import json
import argparse


def assemblePopulation(simtype, infile, popsize):
	if simtype == 'genepop':
		from_genepop = genePop(command_tuple[1])
		virtual_population = (from_genepop[0],
							command_tuple[-2],command_tuple[-1],
							from_genepop[1],from_genepop[2],
							from_genepop[3],from_genepop[4]
							)

	elif simtype == 'simulate':
		virtual_genepool = virtualGenepool(command_tuple[1], command_tuple[2])
		#print virtual_genepool
		number_of_individuals = int(command_tuple[3])
		virtual_population = (virtualPopulation(virtual_genepool, number_of_individuals),
							command_tuple[-2],command_tuple[-1]
							)
	elif simtype == 'simfile':
		input = infile
		with open(input,'r') as o:
			virtual_genepool = json.load(o)
		#print virtual_genepool
		code2loci_dictionary = {}
		for loci_index in range(0,len(virtual_genepool.keys())):
			loci_name = virtual_genepool.keys()[loci_index]
			allele_keys = virtual_genepool[loci_name].keys()
			for allele_code in allele_keys:
				code2loci_dictionary[str(allele_code)] = str(loci_name)
		number_of_individuals = int(popsize)
		# built virtual population from JSON file
		from_virtpop = virtualPopulation(virtual_genepool, number_of_individuals)
		population_count = 1
		population_dictionary = {population_count:[]}
		He_input_dictionary = {}
		sample_dictionary = {}
		for virt_index in range(0,len(from_virtpop.keys())):
			sample_name = from_virtpop.keys()[virt_index]
			genotype_dictionary = from_virtpop[sample_name]
			sample_genotype_list = []
			loci_genotype_list = []
			for allele_index in range(0,len(genotype_dictionary.keys())):
				allele_id = genotype_dictionary.keys()[allele_index]
				allele_list = genotype_dictionary[allele_id]
				for code_index in range(0,len(allele_list)):
					if sample_genotype_list == []:
						sample_genotype_list = [''] * len(allele_list)
					sample_genotype_list[code_index] = str(str(sample_genotype_list[code_index]) +
												str(allele_list[code_index][0]))
			for genotype in sample_genotype_list:
				loci_name = code2loci_dictionary[genotype[:3]]
				if loci_name in He_input_dictionary.keys():
					He_input_dictionary[loci_name].append(genotype)
				else:
					He_input_dictionary[loci_name] = [genotype]
			sample_dictionary[sample_name] = sample_genotype_list
		population_dictionary[population_count] = sample_dictionary
		virtual_population = (sorted(from_virtpop.keys()),
							[str(x) for x in virtual_genepool.keys()], population_count,
							population_dictionary, He_input_dictionary,code2loci_dictionary
							)
	return virtual_population

def virtualGenepool(loci_range, allele_range):
	# creates a virtual pool of loci and alleles to draw from
	loci_min = int(loci_range[0])
	loci_max = int(loci_range[1])
	allele_min = int(allele_range[0])
	allele_max = int(allele_range[1])
	# randomly pick number of loci
	number_of_loci = randint(loci_min,loci_max)
	print 'There are:\n', number_of_loci, ' loci...'
	lociDictionary = {}
	for loci in range(number_of_loci):
		# randomly pick number of alleles
		number_of_alleles = randint(allele_min,allele_max)
		print number_of_alleles, ' alleles for loci',loci,'...'
		max_probability = 100 - (number_of_alleles)
		min_probability = 1
		alleleDone = False
		while not alleleDone == True:
			alleleList = []
			for allele in range(number_of_alleles):
				allele_probability = 1
				if allele + 1 == number_of_alleles:
					allele_probability = max_probability + 1
				elif max_probability > min_probability:
					allele_probability = (allele_probability +
						randint(min_probability,max_probability)
						)
				alleleList.append(float(allele_probability)/float(100))
				max_probability = max_probability - allele_probability + 1
			alleleDone = True
			lociDictionary[loci] = alleleList
	return lociDictionary

def buildGenotype(virtual_genepool):
	#creates
	individual = {0:[], 1:[]}
	for allele_id, alleleList in individual.iteritems():
		for loci, allele_dictionary in virtual_genepool.iteritems():
			maximum = len(allele_dictionary.keys()) - 1
			minimum = 0
			allele_index = randint(minimum, maximum)
			allele_name = allele_dictionary.keys()[allele_index]
			allele_frequency = allele_dictionary[allele_dictionary.keys()[allele_index]]
			probability = float(randint(0, 100)/float(100))
			while probability > allele_frequency:
				allele_index = randint(minimum, maximum)
				allele_name = allele_dictionary.keys()[allele_index]
				allele_frequency = allele_dictionary[allele_dictionary.keys()[allele_index]]		
				probability = float(randint(0, 100)/float(100))
			alleleList.append((allele_name,allele_frequency))
		individual[allele_id] = alleleList
	return individual

def virtualPopulation(virtual_genepool, number_of_individuals):
	sampleList = []
	virtual_population = {}
	for sample in range(number_of_individuals):
		sample_name = 'virtsim_' + str(sample)
		individual = buildGenotype(virtual_genepool)
		sampleList.append(individual)
		virtual_population[sample_name] = individual
	return  virtual_population

def genePop(genepop_file):
	locus_name_list = []
	population_count = 0
	population_dictionary = {}
	virtual_population = []
	He_input_dictionary = {}
	with open(genepop_file,'r+') as file:
		file_contents = file.read().splitlines()
	for line in file_contents:
		line = line.lower()
		if 'title' in line:
			file_description = line.split('\"')[1]
		elif 'pop' not in line:
			if population_count == 0:
				locus_name_list.append(line)
				He_input_dictionary[line] = []
			else:
				sample = line.split(',\t')[0]
				genotypes = line.split(',\t')[1].split('\t')
				for index in range(0,len(genotypes)):
					genotype = genotypes[index]
					He_input_dictionary[locus_name_list[index]].append(genotype)
				sample_pool_id = str(population_count) + '_' + str(sample)
				virtual_population.append(sample_pool_id)
				if population_count in population_dictionary.keys():
					population_dictionary[population_count].append((sample,genotypes))
				else:
					population_dictionary[population_count] = [(sample,genotypes)]
		elif 'pop' in line:
			population_count += 1

	return (virtual_population, locus_name_list, 
			population_count, population_dictionary, 
			He_input_dictionary)

def virtualSeason(virt_pop_tuple):
	population, per_bout, sample_limit, bout_limit, percent_present = virt_pop_tuple
	print len(population),'individuals...'
	print int(percent_present)*100, 'percent are present for each bout...'
	print 'sampling the population...'
	if per_bout != False:
		per_bout = per_bout.split(',')
		season_list = []
		# use this temp list for no resamples
		#temp_list = list(population)
		for n in range(bout_limit):
			print per_bout[n], ' samples taken for bout ', n
			bout_list = []
			temp_list = present(population, percent_present)
			# Warning at this point you must specify more 
			# indivuduals then number of samples per bout
			for m in range(int(per_bout[n])):
				sample = choice(temp_list)
				bout_list.append(sample)
				temp_list.pop(temp_list.index(sample))
			season_list.append(bout_list)
	else:
		print sample_limit, 'samples taken over',bout_limit,'bouts...'
		season_list = []
		# use this temp list for no resamples
		#temp_list = list(population)
		for n in range(bout_limit):
			bout_list = []
			temp_list = present(population, percent_present)
			# Warning at this point you must specify more 
			# indivuduals then number of samples per bout
			for m in range((sample_limit/bout_limit)):
				sample = choice(temp_list)
				bout_list.append(sample)
				temp_list.pop(temp_list.index(sample))
			season_list.append(bout_list)
	return season_list

def present(population, percent_present):
	subset_size = int(len(population))*float(percent_present)
	present_pop = sample(set(population), int(subset_size))
	return present_pop

def resampled(season_sample):
	resampled_list = []
	for x, left in enumerate(season_sample):
		for y, right in enumerate(season_sample):
			common = len(set(left) & set(right))
			if x != y:
				print "bout %s has %s values in common with bout %s"%(x, common, y)
	reference_bout = season_sample[0]
	for bout in season_sample[1:]:
		resampled_list = resampled_list + list(set(reference_bout)&set(bout))
	multiplicity_resampled_dictionary = dict(Counter(resampled_list))
	print len(multiplicity_resampled_dictionary.keys()), 'individuals were resampled...'

	return multiplicity_resampled_dictionary

def geno_cat_calc(He,E1,E2):
	# no dropout
	p0 = (1 - E1)**2
	# single dropout
	p1 = E1*(1 - E1)
	# no false allele
	f0 = (1 - E2)**2
	# single false allele
	f1 = E2*(1 - E2)
	# double false allele
	f2 = E2**2
	# expected frequencies of catigories
	# read as [true genotype]_[observed genotype], ex AA_AA
	AA_AA = ((1 - He)*((p0**2)*(f0**2) + 4*p0*p1*(f0**2) + 
		4*p0*p1*f0*f1 + 4*(p1**2)*(f0**2) + 8*(p1**2)*f0*f1 + 
		4*(p1**2)*(f1**2)) + He*(2*(p1**2)*(f0**2) + 
		4*(p1**2)*f0*f1 + 2*(p1**2)*(f1**2))
		)
	AB_AB = He*((p0**2)*(f0**2))
	AA_AB = ((1 - He)*(4*(p0**2)*f0*f1 + 8*p0*p1*f0*f1 + 
		8*p0*p1*(f1**2)) + He*(4*p0*p1*(f0**2) + 
		8*p0*p1*f0*f1 + 4*p0*p1*(f1**2))
		)
	AA_BB = ((1 - He)*(4*p0*p1*f0*f1 + 4*p0*p1*f0*f2 + 
			8*(p1**2)*f0*f1 + 8*(p1**2)*f0*f2 + 12*(p1**2)*(f1**2)) + 
			He*(2*(p1**2)*(f0**2) + 12*(p1**2)*f0*f1 + 
			14*(p1**2)*(f1**2) + 12*(p1**2)*f0*f2)
			)
	AB_AC = ((1 - He)*(4*(p0**2)*(f1**2)) + He*(4*(p0**2)*f0*f1 + 
			2*(p0**2)*(f1**2))
			)
	AB_CC = ((1 - He)*(2*(p0**2)*f0*f2 + 8*p0*p1*(f1**2) + 
			4*p0*p1*f0*f2) + He*(8*p0*p1*f0*f1 + 
			12*p0*p1*(f1**2) + 8*p0*p1*f0*f2)
			)
	AB_CD = He*(2*(p0**2)*f0*f2 + 2*(p0**2)*(f1**2))
	cat_hz_dictionary = {'AA_AA':AA_AA, 'AB_AB':AB_AB,
						'AA_AB':AA_AB, 'AA_BB':AA_BB,
						'AB_AC':AB_AC, 'AB_CC':AB_CC,
						'AB_CD':AB_CD}
	return cat_hz_dictionary
def He_calc(He_input_dictionary):
	# Calculating He, Var of He, SE, 95% CI of He
	He_locus_dictionary = {}
	for locus, genotypes in He_input_dictionary.iteritems():
		obs_hetero = 0
		total = 0
		allele_dictionary = {}
		for genotype in genotypes:
			if genotype != '000000':
				alleles = genotype[:3], genotype[3:]
				if alleles[0] != alleles[1]:
					obs_hetero += 1
				total += 1
				for allele in alleles:
					if allele in allele_dictionary.keys():
						allele_dictionary[allele] += 1
					else:
						allele_dictionary[allele] = 1
		sum_AF_squared = float(0)
		sum_AF_cubed = float(0)
		total_alleles_counted = float(total*2)
		for allele,count in allele_dictionary.iteritems():
			allele_frequency = float(float(count)/total_alleles_counted)
			sum_AF_squared = sum_AF_squared + float(allele_frequency**2)
			sum_AF_cubed = sum_AF_cubed + float(allele_frequency**3)
			allele_dictionary[allele] = count,allele_frequency

		# PEDANT equations below
		est_He = (total_alleles_counted)*((1-sum_AF_squared)/(total_alleles_counted-1))
		He_variance = ((1/((total_alleles_counted/2)*(total_alleles_counted-1)))*
						(sum_AF_squared-sum_AF_squared**2+2*(total_alleles_counted-2)*
						(sum_AF_cubed-sum_AF_squared**2)))
		He_stddev = float(math.sqrt(He_variance))
		h = He_stddev*1.96
		alpha_95 = est_He-h,est_He+h
		He_locus_dictionary[locus] = (round(est_He,4), round(He_variance,8),
									round(He_stddev,8), round(alpha_95[0],4),
									round(alpha_95[1],4), allele_dictionary)
	return He_locus_dictionary

def cause_FA(sub_genotype, FA_rate_dictionary, locus_name):
	allele_1 = sub_genotype[:3]
	allele_2 = sub_genotype[3:]
	probability = float(randint(0, 100)/float(100))
	FA_rate = float(FA_rate_dictionary[locus_name])
	#print locus_name
	#print probability,FA_rate
	if probability <= FA_rate:
		allele_1 = 'FFF'
	# for allele 2
	probability = float(randint(0, 100)/float(100))
	FA_rate = float(FA_rate_dictionary[locus_name])
	#print probability, FA_rate
	if probability <= FA_rate:
		allele_2 = 'FFF'
	sub_genotype = str(allele_1) + str(allele_2)
	
	return sub_genotype

def cause_ADO(sub_genotype, ADO_rate_dictionary, locus_name):
	allele_1 = sub_genotype[:3]
	allele_2 = sub_genotype[3:]
	if allele_1 != allele_2:
		# for allele 1
		probability = float(randint(0, 100)/float(100))
		ADO_rate = float(ADO_rate_dictionary[locus_name])
		if probability <= ADO_rate:
			allele_1 = allele_2
		# for allele 2
		probability = float(randint(0, 100)/float(100))
		ADO_rate = float(ADO_rate_dictionary[locus_name])
		if probability <= ADO_rate:
			allele_2 = allele_1
		sub_genotype = str(allele_1) + str(allele_2)

	return sub_genotype

def getError(errfile):
	with open(errfile,'r') as e:
		error_data = json.load(e)
		ADO_rate_dictionary = {str(x):float(error_data[x]['ADO']) for x in error_data.keys()}
		FA_rate_dictionary = {str(x):float(error_data[x]['FA']) for x in error_data.keys()}
	return ADO_rate_dictionary, FA_rate_dictionary


###############################################################################################

def main(simtype, infile, errfile, boutlimit,
	samplelimit, perbout, popsize, percentpresent, outfile):
	# assemble the virtual population based on user specifications
	# output object has virt pop, sample_limit, and bout limit
	assemble_it = assemblePopulation(simtype, infile, popsize
		)
	virtual_population = (assemble_it[0], perbout, int(samplelimit), int(boutlimit), percentpresent)
	locus_name_list = assemble_it[1]
	print locus_name_list
	population_dictionary = assemble_it[3]
	# collect virtual samples over a season
	season_sample = virtualSeason(virtual_population)
	# get error rates for each loci from input json file
	ADO_rate_dictionary, FA_rate_dictionary = getError(errfile)

	code2loci_dictionary = assemble_it[5]
	population, perbout, sample_limit, bout_limit, percent_present = virtual_population
	#save_file_name = 'P{0}_SL{1}_BL{2}.csv'.format(unicode(str(len(population)),'utf-8'),
	#											unicode(str(sample_limit),'utf-8'),
	#											unicode(str(bout_limit),'utf-8'))
	save_file_name = outfile
	with open(save_file_name,'w') as o:
		double_locus_name_list = [x for pair in zip(locus_name_list,locus_name_list)
								 for x in pair]
		o.write('sample,' + ','.join(double_locus_name_list) + '\n')
		for bout_index in range(0,len(season_sample)):
			sample_id_list = season_sample[bout_index]
			for sample_id in sample_id_list:
				if sample_id in population_dictionary[1].keys():
					allele_code_list = population_dictionary[1][sample_id]
					new_allele_code_list = []
					for genotype_index in range(0,len(allele_code_list)):
						genotype = allele_code_list[genotype_index]
						locus_name = code2loci_dictionary[genotype[:3]]
						genotype = cause_FA(genotype, FA_rate_dictionary, locus_name)
						genotype = cause_ADO(genotype, ADO_rate_dictionary, locus_name)
						new_allele_code_list.append(genotype)
					code_list_1 = [x[:3] for x in new_allele_code_list]
					code_list_2 = [x[3:] for x in new_allele_code_list]
					double_code_list = [x for pair in zip(code_list_1,code_list_2)
										for x in pair]
					o.write(sample_id + '_' + str(bout_index) + ',' 
							+ ','.join(double_code_list) + '\n')

	# calculate the error probabilities for ADO and FA
	# you must supply heterozygosity, ADO and FA rates for the pop.
	'''
	He_locus_dictionary = He_calc(assemble_it[6])
	#print He_locus_dictionary
	for loci_index in range(0,len(He_locus_dictionary.keys())):
		locus = He_locus_dictionary.keys()[loci_index]
		locus_He = He_locus_dictionary[locus][0]
		ADO = 0.06
		FA = 0.12
		print locus,locus_He
		print '\n'.join([str(x) + ' ' + str(round(y,4)) for x,y in 
						geno_cat_calc(locus_He, ADO, FA).iteritems()])
	'''

	# returns a dictionary of the resampled individuals and their multiplicity
	'''
	resampled_dictionary = resampled(season_sample)
	for sample_id, multiplicity in resampled_dictionary.iteritems():
		print sample_id, 'was resampled',multiplicity,'times...'
	'''


### MAIN ###
if __name__ == '__main__':
# collect arguments from commandline
	parser = argparse.ArgumentParser(description='insert program description')
	parser.add_argument('-t','--simtype', help='Specify which type of simulation type; virtpop, simfile, genepop.', 
		required=True
		)
	parser.add_argument('-i','--infile', help='Specify path to input file; GENEPOP, SIMFILE, or no input.', 
		required=True
		)
	parser.add_argument('-e','--errfile', help='Specify path to input error file; JSON format', 
		required=True
		)
	parser.add_argument('-b','--boutlimit', help='Total number of sampling bouts for a season.', 
		required=True
		)
	parser.add_argument('-s','--samplelimit', help='Total number of samples for a season.', 
		required=True
		)
	parser.add_argument('-l','--perbout', help='Number of samples for each bout.', 
		required=True, default=False
		)
	parser.add_argument('-p','--popsize', help='Specify the size of a virtual population.', 
		required=True
		)
	parser.add_argument('-r','--percentpresent',
		help='Percent of population present at each sampling [default = 100]', 
		required=False, default=100
		)
	parser.add_argument('-o','--outfile', help='Specify the output file.', 
		required=True
		)
	args = vars(parser.parse_args())
	# check to make sure that enough arguments were passed before proceeding
	if len(sys.argv) < 7:
		sys.exit("Missing flags, type \"--help\" for assistance...")
	main(args['simtype'], args['infile'],args['errfile'],
		args['boutlimit'],args['samplelimit'],args['perbout'],args['popsize'],
		args['percentpresent'],args['outfile']
		)








