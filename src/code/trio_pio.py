import sys
import pandas as pd
import itertools
import numpy as np

file = pd.read_csv(sys.argv[1])
model = pd.read_csv(sys.argv[2], sep='\t')
#outfile = sys.argv[3]

upper_bound = model.upper_bound.mean()

#file = file[file['corrected_score'] <= upper_bound]
#samp_tup_list = sorted([sorted(set(x[0] + x[1])) for x in 
#	itertools.combinations(list(zip(file.sample_0,file.sample_1)),2
#	) if len(set(x[0] + x[1])) == 3]
#	)
score_tuples = zip(file.sample_0,file.sample_1,file.corrected_score)
ad_ju_list = []
ju_ta_list = []
ad_ta_list = []
for x in score_tuples:
	if (str(x[0]+x[1]).count('ij') == 1 and (str(x[0]+x[1]).count('m') == 1
		or str(x[0]+x[1]).count('f') == 1)):
		y = tuple(sorted([x[0],x[1]],reverse=True))
		ad_ju_list.append((y[0],y[1],x[2]))
	elif (str(x[0]+x[1]).count('it') == 1 
		and str(x[0]+x[1]).count('ij') == 1):
		y = tuple(sorted([x[0],x[1]]))
		ju_ta_list.append((y[0],y[1],x[2]))
	elif (str(x[0]+x[1]).count('it') == 1 and (str(x[0]+x[1]).count('m') == 1
		or str(x[0]+x[1]).count('f') == 1)):
		y = tuple(sorted([x[0],x[1]],reverse=True))
		ad_ta_list.append((y[0],y[1],x[2]))
ad_ju_df = pd.DataFrame(ad_ju_list, columns=['adults','juveniles','corrected_score'])
ju_ta_df = pd.DataFrame(ju_ta_list, columns=['juveniles','tadpoles','corrected_score'])
ad_ta_df = pd.DataFrame(ad_ta_list, columns=['adults','tadpoles','corrected_score'])
df_list = [ad_ju_df,ju_ta_df,ad_ta_df]
file_names = ['ad_ju_df','ju_ta_df','ad_ta_df']
for df_index in range(0,len(df_list)):
	df = df_list[df_index]
	columns = df.columns.values
	df = df.groupby(columns[0], 
			group_keys=False).apply(lambda x: x.ix[x.corrected_score.idxmin()]
			)

	df = df.groupby(columns[1], 
			group_keys=False).apply(lambda x: x.ix[x.corrected_score.idxmin()]
			)
	df = df.sort(['corrected_score'], ascending=[True]).reset_index(drop=True)
	df.to_csv(file_names[df_index] + '.csv')
	df_list[df_index] = df

trio_df = pd.merge(df_list[0],df_list[1], how='left', on='juveniles')
trio_df = trio_df[trio_df['adults'].isin(df_list[2].adults)]
trio_df = pd.merge(trio_df,df_list[2], how='right', on=['adults','tadpoles'])
trio_df = trio_df.dropna()
trio_df.columns = ['adults','juveniles','ad_ju_score','tadpoles','ju_ta_score','ad_ta_score']
trio_df.to_csv('trio.csv')

'''
ad_ju_score_dictionary = dict((tuple(sorted([x[0],x[1]],reverse=True)),x[2]) for x in score_tuples
	if str(x[0]+x[1]).count('ij') == 1 and (str(x[0]+x[1]).count('m') == 1
	or str(x[0]+x[1]).count('f') == 1))
ju_ta_score_dictionary = dict((tuple(sorted([x[0],x[1]])),x[2]) for x in score_tuples
	if str(x[0]+x[1]).count('it') == 1 and str(x[0]+x[1]).count('ij') ==1)
ad_ta_score_dictionary = dict((tuple(sorted([x[0],x[1]],reverse=True)),x[2]) for x in score_tuples
	if str(x[0]+x[1]).count('it') == 1 and (str(x[0]+x[1]).count('m') == 1
	or str(x[0]+x[1]).count('f') == 1))
'''
'''
with open(outfile,'w') as o:
	for k,v in ad_ju_score_dictionary.iteritems():
		if tup[0] in score_dictionary.keys() and tup[1] in score_dictionary.keys():
			tup_0_score = score_dictionary[tup[0]]
			tup_1_score = score_dictionary[tup[1]]
			#tup = sorted(set(tup[0] + tup[1]))
			tup_string = str(','.join(tup[0]) + ',' + str(tup_0_score)  + ',' 
				+ ','.join(tup[1]) + ',' + str(tup_1_score) + '\n'
				)
			if ('it' in tup_string and 'ij' in tup_string 
				and ('m' in tup_string or 'f' in tup_string)
				):
				o.write(tup_string)
'''