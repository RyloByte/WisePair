###**WisePair**
* * *
####What is it used for?
wisepair is a set of python scripts built to address two issues in individual-based
genetic tracking studies.  
1) Sampling scheme design to maximize resampling a set number of individuals.  
2) Confidently identifying individuals who have been resampled in a data-set.  
* * *
####There are 3 main scripts in the core of WisePair:  
**_beanbag.py_**  
The beanbag.py script was specifically designed to build virtual individual genotypes of a population to be used in simulated sampling.  This design was based on user-supplied criteria such as number of individuals in the population, number of loci, and allelic frequencies.  In addition, this script incorporated genotyping error rates during sampling.  
**_wisepair.py_**  
The second script, wisepair.py, was created to determine the number of re-samples within a specified data set (real or virtual) through allelic pairwise comparisons.  wisepair.py determined the number of re-samples within a virtual data set, determined the number of re-samples within an actual data using specified threshold simulations, estimated the number of errors for re-samples, and determined whether re-samples can be distinguished from non-re-samples.  
**_optimagic.py_**  
The final script, optimagic.py, utilized outputs from both beanbag.py and wisepair.py to develop optimal sampling designs for individual based studies.  beanbag.py and wisepair.py were used to produce a threshold “score” with which we could compare samples to the field data set and subsequent simulations in OPTIMAGIC.py  
* * *
####Usage:
**_beanbag.py_**  
To build a virtual population, beanbag.py requires 2 input files and several flags.  
*****
usage: beanbag.py [-h] -t SIMTYPE -i INFILE -e ERRFILE -b BOUTLIMIT -s  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;SAMPLELIMIT [-l PERBOUT] -p POPSIZE [-r PERCENTPRESENT] -o  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;OUTFILE  

insert program description  

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
*****
