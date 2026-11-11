## require
import sys, pandas, re

## arguments
fin = sys.argv[1] # annotated VCF
fout = sys.argv[2] # annotations extracted into a dataframe


rc = re.compile('transcript\|(.*?)\|protein')

# read file
with open(fin) as flob:
	contents = [ i.strip('\n') for i in flob]

# separate header and contents
headers = [ i for i in contents if i[0]=='#']
contents = [ i for i in contents if i not in headers]
contents = [ i.split('\t') for i in contents]
# extract relevant data into a dataframe
contents = [ [i[0], i[1], i[3]] + [ j for j in i[7].split(';') if 'ANN=' in j] for i in contents]
contents = [ i[:3] + [j] for i in contents for j in i[3].replace('ANN=','').split(',')]
# make sure the transcript names are split into different columns in the next step
for x,i in enumerate(contents):
    mch = re.search( rc, i[-1])
    if bool(mch):
        contents[x][-1] = i[-1].replace( mch.group(1), mch.group(1).replace('|','-'))
contents = [ i[:3] + i[3].split('|') for i in contents]
df = pandas.DataFrame(contents)
# add header
anh = [ i for i in headers if 'ID=ANN' in i][0]
anh = re.search("'(.+)'",anh).group(1).split(' | ')
ogh = headers[-1][1:].split('\t')
ogh = [ ogh[x] for x in [0,1,3]]
df.columns = ogh + anh
df.rename( columns = {'ERRORS / WARNINGS / INFO':'info'},inplace=True)
# separate cds position and length
cpl = df['CDS.pos / CDS.length'].str.split('/',expand=True)
apl = df['AA.pos / AA.length'].str.split('/',expand=True)
df.drop([ 'CDS.pos / CDS.length', 'AA.pos / AA.length', 'cDNA.pos / cDNA.length'],
	axis=1,inplace=True)
df = df.assign( npos=cpl[0], nlen=cpl[1], apos=apl[0], alen=apl[1])
df.replace( '', pandas.NA, inplace=True)
df.to_csv(fout,sep='\t',index=False)
