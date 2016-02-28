* * *
####**Example command to run beanbag.py**
This example builds a virtual population of 300 individuals. It then precedes to sample this population 5 times, totaling 125 samples (25 samples/bout). The -r flag indicates that 50% of the population is present at the start of each bout.  
```python
python beanbag.py -t simfile -i real_data_frequencies.json -e pedant_error_rates.json -b 5 -s 125 -p 300 -r 0.5 -o example_output.csv
```
* * *
####**Example command to run wisepair.py**
This command takes the output from beanbag.py and runs it through wisepair. The output is a scoring matrix for the sample comparisons, a report visualization, and a stats file of the simulations.  
```python
python wisepair.py -s [location of beanbag.py output csv] -l [csv filename with no path]
```
* * *
