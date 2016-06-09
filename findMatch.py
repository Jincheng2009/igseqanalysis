import re 

datapath="C:/Java/data/sarav5/"
CDR = "CDR/"
ref_file = datapath + CDR + "S5.paired.CDR3.count"
ref_data = open(ref_file, "rU")

query_3="KGDL-------------------DV:AAWDSSLSA-------VV"
query_717 = "SPHFDFLGGFPHAS---------DV:ASWDDSLNG-------RV"
query_638 = "GGIDLAPPLGYYYGM--------DV:GTWDSSLSA-------GV"
query = query_3
minMismatch = len(query)
tophit = []

def countMismatch(query, subject, start=0, end=-1):
    if end == -1:
        end = len(query)
    if len(query)!=len(subject):
        return len(query)
    count = 0
    
    if end > len(query):
        end = len(query)
        
    for i in range(start, end):
        if(query[i] != subject[i]):
            count += 1
        #Add penalty for gap
        if query[i] == '-' and subject[i] != '-':
            count += 1
        if query[i] != '-' and subject[i] == '-':
            count += 1
    return count

for line in ref_data:
    tokens = re.split("\t", line)
    seq = tokens[0]
    count = countMismatch(query, seq, 0, 25)
    if count < minMismatch:
        minMismatch = count
        tophit = [seq]
    elif count == minMismatch:
        tophit.append(seq)
        
for item in tophit:
    print item