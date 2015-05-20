from matplotlib import pyplot as plt
import matplotlib.cm as cm
import sys
import pandas as pd
import numpy as np


file = sys.argv[1]
data = pd.read_csv(file)

# separate by min resample limit
# Good RS nd PSC
bouts = list(set(data.number_of_bout))
colors = cm.brg(np.linspace(0, 1, len(bouts)))
legend_series_list = []
for bout in bouts:
    green_data = data.loc[(data.met_resamp_min >= data.min_resampled) &
                          (data.mean_rs_per_mrm >= data.min_times_resampled) & (data.number_of_bout == bout)
                          ]
    # Good RS only
    blue_data = data.loc[(data.met_resamp_min >= data.min_resampled) &
                         (data.mean_rs_per_mrm < data.min_times_resampled) & (data.number_of_bout == bout)
                         ]
    # Good RSC only
    yellow_data = data.loc[(data.met_resamp_min < data.min_resampled) &
                           (data.mean_rs_per_mrm >= data.min_times_resampled) & (data.number_of_bout == bout)
                           ]
    # No criteria met
    red_data = data.loc[(data.met_resamp_min < data.min_resampled) &
                        (data.mean_rs_per_mrm < data.min_times_resampled) & (data.number_of_bout == bout)
                        ]

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

    # red = plt.scatter(rx, ry, c=rc)#, s=rs)
    #blue = plt.scatter(bx, by, c=bc)#, s=bs)
    #yellow = plt.scatter(yx, yy, c=yc)#, s=ys)
    #green = plt.scatter(gx, gy, c=gc)#, s=gs)
    plt.axhline(y=data.min_resampled[0], color="black", ls="--")
    #plt.legend((green, blue, yellow, red),
    #	('RS and RSC Met', 'RS Met only', 'RSC Met only','None Met'),
    #	scatterpoints=1,
    #	loc='upper left'
    #	)

    # Create contours
    contour_data = data.loc[(data.number_of_bout == bout)][['number_of_samples', 'met_resamp_min']]
    h_boundary = contour_data.groupby('number_of_samples')['met_resamp_min'].max().reset_index()
    v_boundary = contour_data.groupby('met_resamp_min')['number_of_samples'].max().reset_index()
    drop_list = []
    for index in range(1, len(h_boundary)):
        found_less = True
        ind = index
        while found_less == True and ind < len(h_boundary) - 1:
            next_row = h_boundary.loc[[ind + 1]]
            current_row = h_boundary.loc[[index]]
            if float(next_row.met_resamp_min) < float(current_row.met_resamp_min):
                drop_list.append(ind + 1)
            if ind == len(h_boundary):
                found_less = False
            ind += 1
    h_boundary.drop(h_boundary.index[drop_list], inplace=True)
    drop_list = []
    for index in range(1, len(v_boundary)):
        found_less = True
        ind = index
        while found_less == True and ind < len(v_boundary) - 1:
            next_row = v_boundary.loc[[ind + 1]]
            current_row = v_boundary.loc[[index]]
            if float(next_row.number_of_samples) < float(current_row.number_of_samples):
                drop_list.append(ind + 1)
            if ind == len(v_boundary):
                found_less = False
            ind += 1
    v_boundary.drop(v_boundary.index[drop_list], inplace=True)
    hori, = plt.plot(h_boundary.number_of_samples, h_boundary.met_resamp_min, color=colors[bouts.index(bout)])
    vert, = plt.plot(v_boundary.number_of_samples, v_boundary.met_resamp_min, color=colors[bouts.index(bout)])
    legend_series_list.append(hori)
plt.legend(legend_series_list,
           bouts,
           loc='upper left',
           ncol=2,
           fontsize=12
           )
# set axis limits
xmin = data.number_of_samples.min() - data.number_of_samples.min() * 0.1
xmax = data.number_of_samples.max() + data.number_of_samples.max() * 0.1
ymin = data.met_resamp_min.min() - 0.5
ymax = data.met_resamp_min.max() + data.met_resamp_min.max() * 0.1

plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.xlabel('total number of samples')
plt.ylabel('# of resampled individuals')
plt.savefig('fig3.pdf')
