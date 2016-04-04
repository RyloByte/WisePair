###**WisePair**
* * *
* * *
###**License**
Copyright (c) 2016 Ryan James McLaughlin
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.
* * *
####What is it used for?
WisePair is a set of python scripts built to address two issues in individual-based
genetic tracking studies.  
1) Sampling scheme design to maximize resampling a set number of individuals.  
2) Confidently identifying individuals who have been resampled in a data-set.  
* * *
####There are 3 main scripts in the core of WisePair:  
**_beanbag.py_**  
The `beanbag.py` script is designed to build virtual individual genotypes of a population to be used in simulated sampling.  This design is based on user-supplied criteria such as number of individuals in the population, number of loci, and allelic frequencies.  In addition, this script incorporates genotyping error rates during sampling.  
**_wisepair.py_**  
The second script, `wisepair.py`, determines the number of re-samples within a specified data set (real or virtual) through allelic pairwise comparisons.  `wisepair.py` determines the number of re-samples within a virtual data set, determines the number of re-samples within an actual data using specified threshold simulations, estimates the number of errors for re-samples, and determines whether re-samples can be distinguished from non-re-samples.  
**_optimagic.py_**  
The final script, `optimagic.py`, utilizes outputs from both `beanbag.py` and `wisepair.py` to develop optimal sampling designs for individual based studies.  `beanbag.py` and `wisepair.py` are used to produce a threshold “score” with which we could compare samples to the field data set and subsequent simulations in `optimagic.py`  
* * *
* * *
####Usage:
####**_beanbag.py_**  
To build a virtual population, `beanbag.py` requires 2 input files and several flags.  
INFILE is a JSON format file which contains loci names, alleles for each loci and allelic frequencies for each allele. ERRFILE is a JSON format file which contains loci names, allelic dropout rates, and false allele rates. Currently there is no function to build these files, so they must be built manually.  
Examples of both files can be found in the **sample_data** directory.  

The required flags are: [-t], [-i], [-e], [-b], [-s], [-p], [-o]  
Optional flags are: [-l], [-r]
An example command is found in the **sample_data** directory, in the example_cmds.txt file.  
* * *
usage: beanbag.py [-h] -t SIMTYPE -i INFILE -e ERRFILE -b BOUTLIMIT -s  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;SAMPLELIMIT [-l PERBOUT] -p POPSIZE [-r PERCENTPRESENT] -o  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;OUTFILE  

optional arguments:  
&nbsp;&nbsp;-h, --help            show this help message and exit  
&nbsp;&nbsp;-t SIMTYPE, --simtype SIMTYPE  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Specify which type of simulation type; currently only  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;simfile  
&nbsp;&nbsp;-i INFILE, --infile INFILE  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Specify path to input frequency file; SIMFILE  
&nbsp;&nbsp;-e ERRFILE, --errfile ERRFILE  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Specify path to input error file; JSON format  
&nbsp;&nbsp;-b BOUTLIMIT, --boutlimit BOUTLIMIT  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Total number of sampling bouts for a season.  
&nbsp;&nbsp;-s SAMPLELIMIT, --samplelimit SAMPLELIMIT  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Total number of samples for a season.  
&nbsp;&nbsp;-l PERBOUT, --perbout PERBOUT  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Number of samples for each bout. [default = False]  
&nbsp;&nbsp;-p POPSIZE, --popsize POPSIZE  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Specify the size of a virtual population.  
&nbsp;&nbsp;-r PERCENTPRESENT, --percentpresent PERCENTPRESENT  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Population present at each sampling between 0-1  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[default = 1]  
&nbsp;&nbsp;-o OUTFILE, --outfile OUTFILE  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Specify the output file.  
* * *
* * *
####**_wisepair.py_**  
`wisepair.py` takes 2 inputs, the directory to the `beanbag.py` output and a list of the output files.  
If real data is being analyzed aswell, those files should be run separatly (and after) from the simulated data so that the model_stats file is already created.  
* * *
usage: wisepair.py [-h] [-s SIMDIR] [-l SIMLIST]  

optional arguments:  
&nbsp;&nbsp;-h, --help            show this help message and exit  
&nbsp;&nbsp;-s SIMDIR, --simdir SIMDIR  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;path to simulation file directory.  
&nbsp;&nbsp;-l SIMLIST, --simlist SIMLIST  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;list of simulation files.  
* * *
* * *
####**_optimagic.py_**  
This script combines requirements of `beanbag.py` and `wisepair.py` along with some ranges that allow for multiple schemes to be tested.  
Real data can be analyzed after models are built for comparison.   
* * *
usage: optimagic.py [-h] [-s SIMDIR] -p POPSIZE [-x PERPOP] -l SAMPRAN -b  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;BOUTRAN -r RESAMPMIN -n NUMMINRAN [-u RUNSIM] [-m ITERLIM]  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[-t SIMTYPE] [-i SIMFILE] [-e ERRFILE]  

optional arguments:  
&nbsp;&nbsp;-h, --help            show this help message and exit  
&nbsp;&nbsp;-s SIMDIR, --simdir SIMDIR  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;path to simulation file directory.  
&nbsp;&nbsp;-p POPSIZE, --popsize POPSIZE  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;population size.  
&nbsp;&nbsp;-x PERPOP, --perpop PERPOP  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;population size present at site.  
&nbsp;&nbsp;-l SAMPRAN, --sampran SAMPRAN  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;range of samples per bout.  
&nbsp;&nbsp;-b BOUTRAN, --boutran BOUTRAN  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;range of bouts.  
&nbsp;&nbsp;-r RESAMPMIN, --resampmin RESAMPMIN  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;minimum individuals resampled over season.  
&nbsp;&nbsp;-n NUMMINRAN, --numminran NUMMINRAN  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;range for minimum number of resamples per resampled  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;individual.  
&nbsp;&nbsp;-u RUNSIM, --runsim RUNSIM  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;use simerator to run simulations.  
&nbsp;&nbsp;-m ITERLIM, --iterlim ITERLIM  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;number of iterations during simulation mode.  
&nbsp;&nbsp;-t SIMTYPE, --simtype SIMTYPE  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;specitfy type of input simulation format: virtpop,  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;simfile, genepop.  
&nbsp;&nbsp;-i SIMFILE, --simfile SIMFILE  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;path to simulation input file.  
&nbsp;&nbsp;-e ERRFILE, --errfile ERRFILE  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Specify path to input error file; JSON format  
* * *
* * *
###Running the simulations to find sampling scheme
* * *
####**Example command to run `beanbag.py`**
This example builds a virtual population of 300 individuals. It then precedes to sample this population 5 times, totaling 125 samples (25 samples/bout). The -r flag indicates that 50% of the population is present at the start of each bout.  
```python
python beanbag.py -t simfile -i real_data_frequencies.json -e pedant_error_rates.json -b 5 -s 125 -p 300 -r 0.5 -o example_output.csv
```
* * *
####**Example command to run `wisepair.py`**
This command takes the output from `beanbag.py` and runs it through wisepair. The output is a scoring matrix for the sample comparisons, a report visualization, and a stats file of the simulations.  
```python
python wisepair.py -s [location of beanbag.py output csv] -l [csv filename with no path]
```
* * *
####**Example command to run `optimagic.py`**
This command runs a optimization run to determine which sampling scheme would work best if 2 individuals need to be resampled at least 3 times from a population of 100 individuals.  
The outputs from `beanbag.py` and `wisepair.py` are sent to the -s path in a directory called **OM_simdir**, and 2 csv files are produced which separate good and bad sampling schemes.  
```python
python optimagic.py -s [path to output] -p 100 -x 0.5 -l 20,25 -b 5,8 -r 2 -n 3 -u True -m 1 -t simfile -i real_data_frequencies.json -e pedant_error_rates.json
```
* * *
* * *
###Analyzing real data with simulations
* * *
####**Example command to run `optimagic.py`**
Run `optimagic.py` with the bounds of the real data sampling scheme at a higher iteration limit to build a model for the real data.  
```python
python optimagic.py -s [path to output] -p 100 -x 0.5 -l 30 -b 5 -r 2 -n 3 -u True -m 10 -t simfile -i real_data_frequencies.json -e pedant_error_rates.json
```
* * *
####**Example command to run `wisepair.py`**
Copy the **model_stats.tsv** from the optimagic output directory to the same location as the real data files.  
Output from this will be identified resampled individuals.  
```python
python wisepair.py -s [location of beanbag.py output csv] -l [csv filename with no path]
```
* * *
