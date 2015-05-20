from matplotlib import pyplot as plt
import sys
import pandas as pd


file = sys.argv[1]
data = pd.read_csv(file)

# separate by min resample limit: g = good, r = bad
green_data = data.loc[(data.met_resamp_min >= data.min_resampled) &
                      (data.mean_rs_per_mrm >= data.min_times_resampled)
                      ]
print green_data.head()

blue_data = data.loc[(data.met_resamp_min >= data.min_resampled) &
                     (data.mean_rs_per_mrm < data.min_times_resampled)
                     ]
print blue_data.head()

yellow_data = data.loc[(data.met_resamp_min < data.min_resampled) &
                       (data.mean_rs_per_mrm >= data.min_times_resampled)
                       ]
print yellow_data.head()

red_data = data.loc[(data.met_resamp_min < data.min_resampled) &
                    (data.mean_rs_per_mrm < data.min_times_resampled)
                    ]
print red_data.head()

gx = green_data.number_of_samples
gy = green_data.met_resamp_min
gs = green_data.number_of_bout
gc = "green"

bx = blue_data.number_of_samples
by = blue_data.met_resamp_min
bs = blue_data.number_of_bout
bc = "blue"

yx = yellow_data.number_of_samples
yy = yellow_data.met_resamp_min
ys = yellow_data.number_of_bout
yc = "yellow"

rx = red_data.number_of_samples
ry = red_data.met_resamp_min
rs = red_data.number_of_bout
rc = "red"

red = plt.scatter(rx, ry, c=rc)  # , s=rs)
blue = plt.scatter(bx, by, c=bc)  # , s=bs)
yellow = plt.scatter(yx, yy, c=yc)  # , s=ys)
green = plt.scatter(gx, gy, c=gc)  # , s=gs)
plt.axhline(y=data.min_resampled[0], color="black", ls="--")
plt.legend((green, blue, yellow, red),
           ('RS and RSC Met', 'RS Met only', 'RSC Met only', 'None Met'),
           scatterpoints=1,
           loc='upper left'
           )
xmin = data.number_of_samples.min() - data.number_of_samples.min() * 0.1
xmax = data.number_of_samples.max() + data.number_of_samples.max() * 0.1
ymin = data.met_resamp_min.min() - 0.5
ymax = data.met_resamp_min.max() + data.met_resamp_min.max() * 0.1

plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.xlabel('total number of samples')
plt.ylabel('# of resampled individuals')
plt.savefig('fig1.pdf')
