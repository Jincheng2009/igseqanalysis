import re 

datapath="C:/Java/data/sarav7/"
CDR = "CDR/"
unpaired = "unpaired/"
ref_file = datapath + CDR + "S3.paired.CDR3.count"

queries = []
queries.append("ELGSI------------------EH:NSRDSRGNH-------VV")
queries.append("DPGGPF-----------------DI:ASWDDSLNS-------WV")
queries.append("DRGYYYM----------------DV:QSYDSSLRG-------SV")
queries.append("KGDL-------------------DV:AAWDSSLSA-------VV")
queries.append("LPLRGYSGYAL------------DY:SSYTSSSTL-------VV")
queries.append("SPHFDFLGGFPHAS---------DV:ASWDDSLNG-------RV")
queries.append("DYSAPYYGLWSWDSPYYYYYTWMDV:QSFDNSLSDS------VT")
queries.append("DYVGMYVSDYYAYGYSPLTYF--DL:YSTDSSGSQ-------GA")
queries.append("DYYGMYSIWGYDLTYRYYYGPYGDY:AAWDDRLSG-------WV")
queries.append("DYSAPYYGLWSWDSPYYYYTWM-DV:QSFDNSLSDS------VT")
queries.append("GGIDLAPPLGYYYGM--------DV:GTWDSSLSA-------GA")


outfile = open("C:/Java/data/sarav7/search/tophit_638_em_S3_VH_full.txt", 'wb')

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

ref_data = open(ref_file, "rU")
minCount = {}
tophit = {}

for query in queries:
    minCount[query] = len(query)

for line in ref_data:
    for query in queries:
        query1 = query #.split(":")[0]
        tokens = re.split("\t", line)
        seq = tokens[0]
        count = countMismatch(query1, seq, 0, 25)
        if count < minCount[query]:
            minCount[query] = count
            tophit[query] = [seq]
        elif count == minCount[query]:
            if query in tophit:
                tophit[query].append(seq)

hitCount = {}
ref_data = open(ref_file, "rU")
allhits = []

# Get the count of hits
for query in tophit:
    allhits.extend(tophit[query])

ref_data = open(ref_file, "rU")    
for line in ref_data:
    tokens = re.split("\t", line)
    seq = tokens[0]
    amount = tokens[1]
    if seq in allhits:
        hitCount[seq] = amount

# Print out the results    
for query in tophit:
    outfile.write("Result for: \n" + query + "\n")
    outfile.write("Tophit(s): \n")
    hits = tophit[query]
    for hit in hits:
        outfile.write(hit + "\t" + str(minCount[query]) + "\t"+ hitCount[hit] + "\n")
    outfile.write("\n")
    
outfile.close()        
