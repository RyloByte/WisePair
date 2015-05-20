from matplotlib import pyplot as plt
import sys
import pandas as pd


file = sys.argv[1]
data = pd.read_csv(file)
data = data.sort(['number_of_samples'], ascending=[True])
x = data.number_of_samples
y1 = data.tot_resamp
y2 = data.met_resamp_min
y3 = data.number_of_bout
blue = plt.plot(x, y1, label='total resampled')
green = plt.plot(x, y2, label='met resamled min')
red = plt.plot(x, y3, label='number of bouts')
plt.legend(loc='upper left')
xmin = data.number_of_samples.min() - data.number_of_samples.min() * 0.1
xmax = data.number_of_samples.max() + data.number_of_samples.max() * 0.1
ymin = data.met_resamp_min.min() - data.met_resamp_min.min() * 0.1
ymax = data.tot_resamp.max() + data.tot_resamp.max() * 0.1

plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.xlabel('total number of samples')
plt.savefig('fig2.pdf')
