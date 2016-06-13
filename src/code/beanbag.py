#! /usr/bin/env python

from random import randint, choice, sample
from collections import Counter
import sys
import math
import json
import argparse


def assemblepopulation(simtype, infile, popsize):
    '''
    Creates a virtual population object with provided
    loci allelic frequency, error rates, and population size.
    '''
    if simtype == 'simfile': # currently only accepts one format for frequencie
        # load input allelic frequenices from JSON file
        input = infile
        with open(input, 'r') as o:
            virtual_genepool = json.load(o)
        # map alleles to their loci
        code2loci_dictionary = {}
        for loci_index in range(0, len(virtual_genepool.keys())):
            loci_name = virtual_genepool.keys()[loci_index]
            allele_keys = virtual_genepool[loci_name].keys()
            for allele_code in allele_keys:
                code2loci_dictionary[str(allele_code)] = str(loci_name)
        # set population limit
        number_of_individuals = int(popsize)
        # built virtual population
        from_virtpop = virtualPopulation(virtual_genepool, number_of_individuals)
        population_count = 1
        population_dictionary = {population_count: []}
        He_input_dictionary = {}
        sample_dictionary = {}
        for virt_index in range(0, len(from_virtpop.keys())):
            sample_name = from_virtpop.keys()[virt_index]
            genotype_dictionary = from_virtpop[sample_name]
            sample_genotype_list = []
            loci_genotype_list = []
            for allele_index in range(0, len(genotype_dictionary.keys())):
                allele_id = genotype_dictionary.keys()[allele_index]
                allele_list = genotype_dictionary[allele_id]
                for code_index in range(0, len(allele_list)):
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
                              population_dictionary, He_input_dictionary, code2loci_dictionary
                              )
        return virtual_population


def virtualPopulation(virtual_genepool, number_of_individuals):
    '''
    Build population from the genepool
    '''

    sampleList = []
    virtual_population = {}
    for sample in range(number_of_individuals):
        sample_name = 'virtsim_' + str(sample)
        # Build individuals genotype
        individual = buildGenotype(virtual_genepool)
        sampleList.append(individual)
        virtual_population[sample_name] = individual
    return virtual_population


def buildGenotype(virtual_genepool):
    '''
    Creates genotype for a single individual from a genepool
    '''
    individual = {0: [], 1: []}
    for allele_id, alleleList in individual.iteritems():
        for loci, allele_dictionary in virtual_genepool.iteritems():
            maximum = len(allele_dictionary.keys()) - 1
            minimum = 0
            allele_index = randint(minimum, maximum)
            allele_name = allele_dictionary.keys()[allele_index]
            allele_frequency = allele_dictionary[allele_dictionary.keys()[allele_index]]
            probability = float(randint(0, 100) / float(100))
            while probability > allele_frequency:
                allele_index = randint(minimum, maximum)
                allele_name = allele_dictionary.keys()[allele_index]
                allele_frequency = allele_dictionary[allele_dictionary.keys()[allele_index]]
                probability = float(randint(0, 100) / float(100))
            alleleList.append((allele_name, allele_frequency))
        individual[allele_id] = alleleList
    return individual

def virtualSeason(virt_pop_tuple):
    '''
    Run a sampling season on the virtual population
    '''
    population, per_bout, sample_limit, bout_limit, percent_present = virt_pop_tuple
    print len(population), 'individuals...'
    print float(percent_present) * 100, 'percent are present for each bout...'
    print 'sampling the population...'
    if per_bout != 'False':
        per_bout = per_bout.split(',')
        season_list = []
        # use this temp list for no resamples
        # temp_list = list(population)
        for n in range(bout_limit):
            print per_bout[n], ' samples taken for bout ', n
            bout_list = []
            temp_list = present(population, percent_present)
            # Warning at this point you must specify more 
            # indivuduals then number of samples per bout
            for m in range(int(per_bout[n])):
                sample = choice(temp_list)
                bout_list.append(sample)
                # changed this so that it assumes replacement, for pop < sample number per bout
                # temp_list.pop(temp_list.index(sample))
            season_list.append(bout_list)
    else:
        print sample_limit, 'samples taken over', bout_limit, 'bouts...'
        season_list = []
        # use this temp list for no resamples
        # temp_list = list(population)
        for n in range(bout_limit):
            bout_list = []
            temp_list = present(population, percent_present)
            # Warning at this point you must specify more 
            # indivuduals then number of samples per bout
            for m in range((sample_limit / bout_limit)):
                sample = choice(temp_list)
                bout_list.append(sample)
                temp_list.pop(temp_list.index(sample))
            season_list.append(bout_list)
    return season_list


def present(population, percent_present):
    '''
    Build sub-population of individuals present at the time of sampling
    '''
    subset_size = int(len(population)) * float(percent_present)
    present_pop = sample(set(population), int(subset_size))
    return present_pop


def resampled(season_sample):
    '''
    
    '''
    resampled_list = []
    for x, left in enumerate(season_sample):
        for y, right in enumerate(season_sample):
            common = len(set(left) & set(right))
            if x != y:
                print "bout %s has %s values in common with bout %s" % (x, common, y)
    reference_bout = season_sample[0]
    for bout in season_sample[1:]:
        resampled_list = resampled_list + list(set(reference_bout) & set(bout))
    multiplicity_resampled_dictionary = dict(Counter(resampled_list))
    print len(multiplicity_resampled_dictionary.keys()), 'individuals were resampled...'

    return multiplicity_resampled_dictionary


def cause_FA(sub_genotype, FA_rate, locus_name):
    allele_1 = sub_genotype[:3]
    allele_2 = sub_genotype[3:]
    probability = float(randint(0, 100) / float(100))
    # print locus_name
    # print probability,FA_rate
    if probability <= FA_rate:
        allele_1 = 'FFF'
    # for allele 2
    probability = float(randint(0, 100) / float(100))
    # print probability, FA_rate
    if probability <= FA_rate:
        allele_2 = 'FFF'
    sub_genotype = str(allele_1) + str(allele_2)

    return sub_genotype


def cause_ADO(sub_genotype, ADO_rate, locus_name):
    allele_1 = sub_genotype[:3]
    allele_2 = sub_genotype[3:]
    rep_1 = 1
    rep_2 = 1
    if allele_1 != allele_2:
        # for allele 1
        probability = float(randint(0, 100) / float(100))
        if probability <= ADO_rate:
            allele_1 = allele_2
        # for allele 2
        probability = float(randint(0, 100) / float(100))
        if probability <= ADO_rate:
            allele_2 = allele_1
        sub_genotype = str(allele_1) + str(allele_2)

    return sub_genotype


def getError(errfile):
    with open(errfile, 'r') as e:
        error_data = json.load(e)
        ADO_rate_dictionary = {str(x): float(error_data[x]['ADO']) for x in error_data.keys()}
        FA_rate_dictionary = {str(x): float(error_data[x]['FA']) for x in error_data.keys()}
    return ADO_rate_dictionary, FA_rate_dictionary


###############################################################################################

def main(simtype, infile, errfile, boutlimit,
         samplelimit, perbout, popsize, percentpresent, outfile):
    # assemble the virtual population based on user specifications
    # output object has virt pop, sample_limit, and bout limit
    assemble_it = assemblepopulation(simtype, infile, popsize)
    virtual_population = (assemble_it[0], str(perbout), int(samplelimit), int(boutlimit), percentpresent)
    locus_name_list = assemble_it[1]
    print locus_name_list
    population_dictionary = assemble_it[3]
    # collect virtual samples over a season
    season_sample = virtualSeason(virtual_population)
    # get error rates for each loci from input json file
    ADO_rate_dictionary, FA_rate_dictionary = getError(errfile)

    code2loci_dictionary = assemble_it[5]
    population, perbout, sample_limit, bout_limit, percent_present = virtual_population
    # save_file_name = 'P{0}_SL{1}_BL{2}.csv'.format(unicode(str(len(population)),'utf-8'),
    #                                            unicode(str(sample_limit),'utf-8'),
    #                                            unicode(str(bout_limit),'utf-8'))
    '''
    save_file_name = outfile
    with open(save_file_name, 'w') as o:
        double_locus_name_list = [x for pair in zip(locus_name_list, locus_name_list)
                                  for x in pair]
        o.write('sample,' + ','.join(double_locus_name_list) + '\n')
        for bout_index in range(0, len(season_sample)):
            sample_id_list = season_sample[bout_index]
            for sample_id in sample_id_list:
                if sample_id in population_dictionary[1].keys():
                    allele_code_list = population_dictionary[1][sample_id]
                    new_allele_code_list = []
                    for genotype_index in range(0, len(allele_code_list)):
                        genotype = allele_code_list[genotype_index]
                        locus_name = code2loci_dictionary[genotype[:3]]
                        genotype = cause_FA(genotype, FA_rate_dictionary, locus_name)
                        genotype = cause_ADO(genotype, ADO_rate_dictionary, locus_name)
                        new_allele_code_list.append(genotype)
                    code_list_1 = [x[:3] for x in new_allele_code_list]
                    code_list_2 = [x[3:] for x in new_allele_code_list]
                    double_code_list = [x for pair in zip(code_list_1, code_list_2)
                                        for x in pair]
                    o.write(sample_id + '_' + str(bout_index) + ','
                            + ','.join(double_code_list) + '\n')
    '''
    # EXPERIMENTAL: trying to incorporate a repeat function to deal with high ADO error rates                
    save_file_name = outfile
    with open(save_file_name, 'w') as o:
        double_locus_name_list = [x for pair in zip(locus_name_list, locus_name_list)
                                  for x in pair]
        o.write('sample,' + ','.join(double_locus_name_list) + '\n')
        for bout_index in range(0, len(season_sample)):
            sample_id_list = season_sample[bout_index]
            for sample_id in sample_id_list:
                if sample_id in population_dictionary[1].keys():
                    allele_code_list = population_dictionary[1][sample_id]
                    new_allele_code_list = []
                    for genotype_index in range(0, len(allele_code_list)):
                        genotype = allele_code_list[genotype_index]
                        locus_name = code2loci_dictionary[genotype[:3]]
                        FA_rate = float(FA_rate_dictionary[locus_name])
                        ADO_rate = float(ADO_rate_dictionary[locus_name])
                        genotype = cause_FA(genotype, FA_rate, locus_name)
                        # EXPERIMENTAL: dealing with high ADO error rate
                        # First I try will be reducing error to under 10%
                        if ADO_rate > 0.1:
                            print 'error rate to high, you need repeats...'
                            repeat = 1
                            stop = False
                            while stop is False:
                                new_genotype = cause_ADO(genotype, ADO_rate, locus_name)
                                if new_genotype != genotype:
                                    new_genotype = cause_ADO(genotype, ADO_rate, locus_name)
                                    repeat += 1
                                else:
                                    stop = True
                        else:
                            genotype = cause_ADO(genotype, ADO_rate, locus_name)
                        print locus_name, genotype, repeat
                        new_allele_code_list.append(genotype)
                    code_list_1 = [x[:3] for x in new_allele_code_list]
                    code_list_2 = [x[3:] for x in new_allele_code_list]
                    double_code_list = [x for pair in zip(code_list_1, code_list_2)
                                        for x in pair]
                    o.write(sample_id + '_' + str(bout_index) + ','
                            + ','.join(double_code_list) + '\n')


### MAIN ###
if __name__ == '__main__':
    # collect arguments from commandline
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-t', '--simtype', help='Specify which type of simulation type; currently only simfile',
                        required=True
                        )
    parser.add_argument('-i', '--infile', help='Specify path to input frequency file; SIMFILE',
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
    parser.add_argument('-l', '--perbout', help='Number of samples for each bout. [default = False]',
                        required=False, default='False'
                        )
    parser.add_argument('-p', '--popsize', help='Specify the size of a virtual population.',
                        required=True
                        )
    parser.add_argument('-r', '--percentpresent',
                        help='Population present at each sampling between 0-1 [default = 1]',
                        required=False, default=1
                        )
    parser.add_argument('-o', '--outfile', help='Specify the output file.',
                        required=True
                        )
    args = vars(parser.parse_args())
    # check to make sure that enough arguments were passed before proceeding
    if len(sys.argv) < 7:
        sys.exit("Missing flags, type \"--help\" for assistance...")
    main(args['simtype'], args['infile'], args['errfile'],
         args['boutlimit'], args['samplelimit'], args['perbout'], args['popsize'],
         args['percentpresent'], args['outfile']
         )
