* * *
####**Example command to run beanbag.py**
This example builds a virtual population of 300 individuals. It then precedes to sample this population 5 times, totaling 125 samples (25 samples/bout). The -r flag indicates that 50% of the population is present at the start of each bout.  
```python
python beanbag.py -t simfile -i real_data_frequencies.json -e pedant_error_rates.json -b 5 -s 125 -p 300 -r 0.5 -o example_output.csv
```
* * *

