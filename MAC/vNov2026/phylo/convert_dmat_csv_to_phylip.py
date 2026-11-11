### convert a CSV-format distance matrix to Phylip 
import sys, pandas

fin = sys.argv[1]
fout = sys.argv[2]

df = pandas.read_csv( fin, index_col=0)
ns = df.shape[0]

for i in range(ns-1):
	for j in range(i+1,ns):
		df.iloc[j,i] = df.iloc[i,j]

for i in range(ns):
	df.iloc[i,i] = 0

df = df.round(decimals=5)

with open(fout,'w') as flob:
	flob.write( str(ns)+'\n')
	df.to_csv( flob, sep='\t', header=False,float_format='{:6.5f}'.format)
