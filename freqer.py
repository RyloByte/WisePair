import pandas as pd
import sys
import json

file = sys.argv[1]
outfile = sys.argv[2]
df = pd.read_csv(file)
print df.columns
for column in df.columns:
	if column.find('.') != -1:
		df[column.split('.')[0]] = df[column.split('.')[0]] + df[column]
		del df[column]
del df['sample']
loci_dict = {}
for column in df.columns:
	count_up = dict(df[column].value_counts())
	count_up.pop(0,None)
	sum_up = sum(count_up.values())
	for alle,count in count_up.iteritems():
		freq = round(float(count)/float(sum_up),3)
		count_up[alle] = str(freq)
	loci_dict[column] = count_up
with open(outfile,'w') as o:
	json.dump(loci_dict,o)
