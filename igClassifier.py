import pandas as pd
import re
import cx_Oracle
import xml.etree.ElementTree as ET
from pandas import DataFrame
import numpy as np


class Sequence:
# Common base class for sequence with annotation
    def __init__(self, sequence, name):
        self.name = name
        self.sequence = sequence
        self.annotations = []
        self.CDR= []
        
    def addRangeAnnotation(self, ra):
        self.annotations.append(ra)  
        if(re.match('^(Light|Heavy) chain CDR',ra.getType())):
            self.CDR.append(ra)

    def getAllAnnotations(self):
        return self.annotations
    
    def getAnnotation(self, position):
        local=[]
        for ra in self.annotations:
            if(position >= ra.getRange()[0] and position<=ra.getRange[1]):
                local.append(ra)
        return local
    
    # return 2 if out of range
    # return 0 if in CDR
    # return 1 if not in CDR
    def isCDR(self, position):
        #print position
        if(position <0 or position >= len(self.sequence)):
            return 2
        else:
            for ra in self.CDR:
                #print (ra.start, ra.end)
                if(position >= ra.start and position<=ra.end):
                    return 0
        return 1
            
class RangeAnnotation:
# Common base class for annotation
    def __init__(self, ratype, value, start, end):
        self.ratype = ratype
        self.value = value
        self.start = start
        self.end = end
        
    def getRange(self):
        rarange = [self.start, self.end]
        return rarange
    
    def getType(self):
        return self.ratype
        
    def getValue(self):
        return self.value
                
# Make connection to qa database
qa_con = cx_Oracle.connect("seqreq/reqseq@uk2as122:1541/qa.cambridgeantibody.com")
cursor = qa_con.cursor()

# Fetch the data from QA database
fetchGermline = """SELECT gl.seq, seq.seqdat, ba.xmldata, gl.name
              FROM germline gl,seq, blaze_annot ba 
              WHERE gl.species='hs' and gl.name LIKE 'IG%' 
              AND gl.seq = seq.id 
              AND ba.seq = seq.id """
              
cursor.execute(fetchGermline)

germlineDict = {}
for row in cursor: 
    blazexml = str(row[2])
    name = row[3]
    root = ET.fromstring(blazexml)
    sequence = root[0].find('sequence').text
    germline = Sequence(sequence, name)
    print "Caching annotation of germline :" + name  
    for ra in root[0][1].findall('ra'):
        ratype = ra.find('t').text
        value = ra.find('v').text
        start = ra.find('s').text
        end = ra.find('e').text
        germline.addRangeAnnotation(RangeAnnotation(ratype, value, int(start), int(end)))
    germlineDict[name] = germline

# Close connection
cursor.close()
qa_con.close()


# Test
for ra in germlineDict['IGHJ1*01'].getAllAnnotations():
    print ra.getType() + "\t" + ra.getValue() + "\t" + str(ra.getRange()[0]) + "\t" + str(ra.getRange()[1])
    
######################################
## Custom import data section
######################################
##########Andrew's data
datapath="C:/java/data/Andrew/mutation/"
igM_0="0-IgM_S13_L001_R1_001.vh.mutation.csv"
igG_0="0-IgG_S21_L001_R1_001.vh.mutation.csv"
# igA_24="24-IgA_S23.mutation.csv"
  
df = pd.read_csv(datapath + igM_0, header=None)
df['isotype'] = 'IgM'
 
temp = pd.read_csv(datapath + igG_0, header=None)
temp['isotype'] = 'IgG'
df = pd.concat([df, temp])
 
# temp = pd.read_csv(datapath + igA_24, header=None)
# temp['isotype'] = 'IgA'
# df = pd.concat([df, temp])
 
df.columns = ['id', 'germline', 'query_pos', 'ref_pos', 'query_base', 'ref_base', 'b3', 'b2', 'b1', 'a1', 'a2', 'a3', 'mismatch', 'length', 'gap', 'phred', 'refstart','refend','isotype']

########Steven's data
datapath="C:/java/data/steven/mutation/"
sample1="SRR1383461_2.vh.mutation.csv"
# sample1="SRR1383461.mutation.csv"
# sample2="SRR1383463.mutation.csv"
# sample3="SRR1383466.mutation.csv"
#   
df = pd.read_csv(datapath + sample1, header=None)
df['sample'] = 'SRR1383461'
df = df.iloc[0:1000000,]
#  
# temp = pd.read_csv(datapath + sample2, header=None)
# temp['sample'] = 'SRR1383463'
# df = pd.concat([df, temp])
#  
# temp = pd.read_csv(datapath + sample3, header=None)
# temp['sample'] = 'SRR1383466'
# df = pd.concat([df, temp])
#  
df.columns = ['id', 'germline', 'query_pos', 'ref_pos', 'query_base', 'ref_base', 'b3', 'b2', 'b1', 'a1', 'a2', 'a3', 'mismatch', 'length', 'gap', 'phred', 'sample']
# df['id'] = df['sample'].astype(str) + "." +df['id'].astype(str)
isotype=pd.read_csv(datapath+"isotype.csv")
isotype = isotype[['id','isotype']]

df['sid'] = df['sample'] + "." + df['id'].map(str)
df = df.drop(['id'],axis=1)
df = df.rename(columns={'sid':'id'})
df = pd.merge(df, isotype, on='id')
######################################
######################################
df = df.drop(['b3','b2','b1','a1','a2','a3'], axis=1)
df['germline'] = df['germline'].str.strip()

###########Test
mut_df = df.drop_duplicates(['id','germline','mismatch','length'])
mut_df1 = mut_df.groupby(['id', 'isotype']).agg({
                                     'length': np.max
                                     })

mut_df1 = mut_df1.reset_index()
mut_df2 = mut_df[['id','germline','length','mismatch']]

mut_df3 = pd.merge(mut_df1, mut_df2, on=['id','length'])

mut_df3['misfrac'] = mut_df3['mismatch']/mut_df3['length']
mut_df3.groupby(['isotype'])['misfrac'].mean()

########################################End test   
germlineTab={}   
for key in germlineDict:
    coordinate ={}
    gene = germlineDict[key]
    for cdr in gene.CDR:
        if "CDR1" in cdr.getType():
            coordinate['cdr1_start']=cdr.getRange()[0]
            coordinate['cdr1_end']=cdr.getRange()[1]
        elif "CDR2" in cdr.getType():
            coordinate['cdr2_start']=cdr.getRange()[0]
            coordinate['cdr2_end']=cdr.getRange()[1]
    coordinate['ref_length'] = len(gene.sequence)
    germlineTab[key] = coordinate      
      
germline_df = DataFrame.from_dict(germlineTab, orient="index")
germline_df['germline'] = list(germline_df.index)

df2 = pd.merge(df, germline_df, on='germline')

def inVseg(x):
    if x['ref_pos'] >= x['ref_length'] or x['ref_pos'] <0:
        return 0
    else:
        return 1
    
def inCDR1(x):
    if x['ref_pos'] >= x['cdr1_start'] and x['ref_pos'] <= x['cdr1_end']:
        return 1
    else:
        return 0
    
def inCDR2(x):
    if x['ref_pos'] >= x['cdr2_start'] and x['ref_pos'] <= x['cdr2_end']:
        return 1
    else:
        return 0
    
columnTemp = df2.apply(inVseg, axis=1)
df2['inVseg'] = columnTemp
columnTemp = df2.apply(inCDR1, axis=1)
df2['inCDR1'] = columnTemp
columnTemp = df2.apply(inCDR2, axis=1)
df2['inCDR2'] = columnTemp



# Temporary save the data
# df2.to_csv(datapath+"0-IgM_S13.mutation2.csv")
df3 = df2.drop_duplicates(['id','query_pos'])

# Subset only the v-gene mutations
v_df = df3[df3['germline'].str.startswith('IGHV')]
j_df = df3[df3['germline'].str.startswith('IGHJ')]

grouped = v_df.groupby(['id', 'germline'])

v_mut_df = grouped.agg({
                      'inVseg': np.sum,
                      'inCDR1': np.sum,
                      'inCDR2': np.sum,
                      'length': np.max
                      })

v_mut_df = v_mut_df.reset_index()

j_mut_df = j_df.groupby(['id', 'germline']).agg({
                                                 'inVseg': np.size,
                                                 'length': np.max
                                                 })

j_mut_df = j_mut_df.reset_index()

mut_df = pd.merge(v_mut_df, j_mut_df, on='id')

mut_df.columns = ['id', 'vgene', 'inCDR1', 'inCDR2', 'vlength', 'vcount', 'jgene', 'jlength', 'jcount']
mut_df['vrate'] = mut_df['vcount'] / mut_df['vlength']
mut_df['jrate'] = mut_df['jcount'] / mut_df['jlength']

phred_df = df.groupby(['id','query_pos'])['phred'].max()
phred_df = phred_df.reset_index()
err_df = phred_df.groupby('id')['phred'].mean()
err_df = err_df.reset_index()

mut_df = pd.merge(mut_df, err_df, on = 'id')


mut_df.to_csv(datapath+"Ig.mutation.csv")
