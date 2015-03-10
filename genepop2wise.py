import sys

file = sys.argv[1]
outfile = file.rsplit('.')[0] + '.wise'

with open(file,'r') as o:
	data = [x.split('\r\n')[0] for x in o.readlines()]
	half_loci = [x for x in data[1:data.index('pop')]]
	loci_names = [None]*(len(half_loci)+len(half_loci))
	loci_names[::2] = half_loci
	loci_names[1::2] = half_loci
	population = [(x.split(' , ',1)[0],x.split(' , ',1)[1]) for x in data[data.index('pop') + 1:]]
with open(outfile,'w') as o:
	o.write('sample,'  + ','.join(loci_names) + '\n')
	for sample_id,geno in population:
		out_list = []
		out_list.append(sample_id)
		geno_split = geno.split(' ')
		for alle in geno_split:
			a1 = str(alle[:3])
			a2 = str(alle[3:])
			out_list.append(a1)
			out_list.append(a2)
		o.write(','.join(out_list) + '\n')
