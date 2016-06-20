import sys

complement={}
complement['A'] = 'T'
complement['T'] = 'A'
complement['G'] = 'C'
complement['C'] = 'G'
complement['-'] = '-'
complement['N'] = 'N'
complement['NA'] = 'NA'

bases = ['A','T','C','G']

codontable = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'O', 'TAG':'B',
'TGC':'C', 'TGT':'C', 'TGA':'J', 'TGG':'W',
}


class Alignment:
    
    def __init__(self, name, length):
        self.name=name
        self.length = length
        self.query = None
        self.vseq = None
        self.jseq = None
        
    def setQuery(self, seq):
        self.query=seq
        
    def hasQuery(self,seq):
        if seq in self.alignDict:
            return True
        return False
        
    def add_v_alignment(self, reference, start):
        self.vseq = reference
        pad = "".join('-' for i in range(start))
        self.vseq.addPrefix(pad)
    
    def add_j_alignment(self, reference, start):
        self.jseq = reference
        pad = "".join('-' for i in range(start))
        self.jseq.addPrefix(pad)
     
    def getCurrentQuery(self):
        return self.query   
        
    def printOut(self):
        if self.query.getTranslation() is not None:
            sys.stdout.write(self.query.getTranslation() + "\n")
        sys.stdout.write(self.query.getSequence() + "\n")
        if self.vseq is not None:
            sys.stdout.write(self.vseq.getSequence() + "\n")
        if self.jseq is not None:
            sys.stdout.write(self.jseq.getSequence() + "\n")
    
    def get_v_sequence(self):
        return self.vseq
    
    def get_j_sequence(self):
        return self.jseq
                
    def getMutations(self):
        mutations = []
        query = self.query
        query_seq = query.getSequence()
        refs = []
        if self.vseq is not None:
            refs.append(self.vseq)
        if self.jseq is not None:
            refs.append(self.jseq)
        for seq in refs:
            rstart = seq.getRange()[0]
            ref_seq = seq.getSequence()
            germline = seq.getName()
            for j in range(len(ref_seq)):
                # when a mismatch is detected
                if ref_seq[j] in bases and query_seq[j] in bases:
                    ref_pos = rstart + seq.getUngappedPosition(j)
                    # Get the two nucleotides before and after
                    query_b1="-"
                    query_b2="-"
                    query_a1="-"
                    query_a2="-"
                    if j >=2:
                        query_b1 = query_seq[j-1]
                        query_b2 = query_seq[j-2]
                    if j < len(query_seq) - 2:
                        query_a1 = query_seq[j+1]
                        query_a2 = query_seq[j+2]  
                    # Get the two nucleotides before and after
                    ref_b1=query_b1
                    ref_b2=query_b2
                    ref_a1=query_a1
                    ref_a2=query_a2
                    if j >=2:
                        if ref_seq[j-1] != ".":
                            ref_b1 = ref_seq[j-1]
                        if ref_seq[j-2] != ".":
                            ref_b2 = ref_seq[j-2]
                    if j < len(query_seq) - 2:
                        if ref_seq[j+1] != ".":
                            ref_a1 = ref_seq[j+1]
                        if ref_seq[j+2] != ".":
                            ref_a2 = ref_seq[j+2]     
                    ref_base = ref_seq[j]
                    query_base = query_seq[j]
                    query_position = query.getRange()[0] + query.getUngappedPosition(j)
                    if query.getTranslation() is not None:
                        prot = query.getTranslation()
                        query_codon=""
                        ref_codon=""
                        if prot[j] != " ":
                            query_codon =  query_b1 + query_base + query_a1
                            ref_codon = ref_b1 + ref_base + ref_a1
                        elif (j-1)>=0 and prot[j-1] != " ":
                            query_codon = query_b2 + query_b1 + query_base
                            ref_codon = ref_b2 + ref_b1 + ref_base
                        elif (j+1)<len(ref_seq) and prot[j+1] != " ":
                            query_codon = query_base + query_a1 + query_a2
                            ref_codon = ref_base + ref_a1 + ref_a2
                        if query_codon in codontable:
                            aa = codontable[query_codon]
                        else:
                            aa = "-"
                        if ref_codon in codontable:
                            ref_aa =  codontable[ref_codon]
                        else:
                            ref_aa = "-"
                    if not query.isStrandPlus():
                        query_position = self.length - query_position + 1
                    record=[self.name, germline, query_position, ref_pos, query_base, ref_base, aa, ref_aa]
                    record.extend([query_b2, query_b1, query_a1, query_a2, ref_b2, ref_b1, ref_a1, ref_a2])
                    mutations.append(record)           
        return mutations                

class Sequence:
    
    def __init__(self, seq, start, end, name = None, strandPlus=True):
        self.seq = seq
        self.start = int(start)
        self.end = int(end)
        self.strandPlus=strandPlus
        self.name = name
        # if gap=[40,3], it means that there are 3 gaps before position 40 (p<=40)
        self.gap=[0,0]
        self.gap[0]=0
        self.prot=""
        
    def getRange(self):
        return [int(self.start), int(self.end)]
    
    def getSequence(self):
        return self.seq
    
    def getName(self):
        return self.name
    
    def addSequence(self, seq, start, end):
        if end > self.end:
            self.end = end
        self.seq = self.seq + seq
    
    def addPrefix(self, seq):
        self.seq = seq + self.seq        
    
    def getUngappedPosition(self, p):
        ngap = 0
        for i in range(0, p+1):
            if self.seq[i] == "-":
                ngap = ngap + 1          
        ungap_pos = p - ngap
        return ungap_pos

    def addTranslation(self, prot):
        self.prot=self.prot + prot

    def getTranslation(self):
        return self.prot
    
    def isStrandPlus(self):
        return self.strandPlus
       