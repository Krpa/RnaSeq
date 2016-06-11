import sys

CHUNK_LENGTH = 60
T_ID = "transcript_id"
LOOK_S = False
ALT = False

class Exon:
    """
    Represents an exon parsed from a GTF file.
    https://genome.ucsc.edu/FAQ/FAQformat.html#format4
    """
    def __init__(self, seqname, source, feature, start, end, score, strand, frame, attribute):
        self.seqname = seqname.split(" ")[0]
        self.source = source
        self.feature = feature
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.frame = frame
        self.attribute = attribute

        for a in attribute.split(";"):
            key,value = a.strip().split(" ")
            if key == T_ID:
                self.trans_id = value.strip('"')
                break

        if self.trans_id is None:
            print "No transcript id!\n" + self

    def __str__(self):
        return "Exon: \n" + " ".join([self.seqname, self.source, self.feature, \
            str(self.start), str(self.end), self.score, self.strand, self.frame, self.attribute])

def complement(base):
    if(base == 'A'):
        return 'T'
    if(base == 'T'):
        return 'A'
    if(base == 'G'):
        return 'C'
    if(base == 'C'):
        return 'G'
    return 'X'

def complementString(string):
    s = ""
    for c in string:
        s = s + complement(c)
    return s

def parse(lines):
    """
    Parser for a GTF file. Extracts all exons and returns 2 dicts.
    First dict returned maps transcript id to a list of exons.
    Second dict returned maps transcript id to sequence id.
    """
    exons = {}
    transToSeq = {}
    for line in lines:
        fields = line.split("\t")
        if len(fields) < 9 or fields[2] != "exon":
            continue
        exon = Exon(*fields)
        if exon.trans_id is None:
            continue
        if exon.trans_id not in exons:
            exons[exon.trans_id] = []
        exons[exon.trans_id].append(exon)
        if exon.trans_id not in transToSeq:
            transToSeq[exon.trans_id] = (exon.seqname, exon.strand)
    return exons, transToSeq

def makeRegions(tid_exons):
    """
    Resolves lists of exons in such way that it combines overlapping
    exons into regions.
    Returns dict that maps transcript id to list of regions.
    """
    tid_regions = {}
    for tid in tid_exons:
        exons = sorted(tid_exons[tid], key=lambda x: x.start)
        start = exons[0].start
        end = exons[0].end
        strand = exons[0].strand
        regions = []
        for idx, exon in enumerate(exons):
            if exon.start <= end:
                end = max(end, exon.end)
            else:
                regions.append((start, end))
                start = exon.start
                end = exon.end
        if len(regions) == 0 or regions[-1] != (start, end):
            regions.append((start, end))
        tid_regions[tid] = regions
    return tid_regions

def altSplicing(regions):
    """
    Generator for alternative splicing of given exons.
    Yields a list of regions and a suffix representing the bitmask
    used for splicing.
    """
    M = len(regions)
    N = (1 << M) - 1
    for i in range(1, N):
        suff = ""
        nextRegions = []
        for j in range(0, M):
            if i & (1 << j) != 0:
                suff += '1'
                nextRegions.append(regions[j])
            else:
                suff += '0'
        yield (nextRegions, suff)

plus = 0
minus = 0
duljina = 0
def makeTranscript(seq, regions, strand):
    """
    Extracts all regions from a sequence.
    """
    global plus
    global minus
    global duljina
    t = ""
    if strand == '+' or not LOOK_S:
        plus += 1
        for start, end in regions:
            t = t + seq[start-1:end]
    else:
        minus += 1
        for start, end in reversed(regions):
            t = t + complementString(seq[start-1:end][::-1])
    duljina += len(t)
    return t

def chunks(string):
    """
    Generator for string slicing.
    """
    return (string[0+i:CHUNK_LENGTH+i] for i in range(0, len(string), CHUNK_LENGTH))

def writeTranscript(out_fasta, header, t):
    """
    Writes transcript to fasta file.
    """
    out_fasta.write(header + '\n')
    for chunk in chunks(t):
        out_fasta.write(chunk + '\n')

def solveSeq(out_fasta, header, seq, toSolve, tid_regions, transToSeq):
    """
    Matches transcripts to given sequence, extracts them and
    writes them to a file.
    """
    name = header[1::].split(" ")[0]
    transcripts = [ind for ind in toSolve if name == transToSeq[ind][0]]
    for t in transcripts:
        tm = makeTranscript(seq, tid_regions[t], transToSeq[t][1])
        writeTranscript(out_fasta, ">" + t, tm)
        if ALT:
            for regs, suff in altSplicing(tid_regions[t]):
                tm = makeTranscript(seq, regs, transToSeq[t][1])
                writeTranscript(out_fasta, ">" + t + "_" + suff + " A", tm)
    return [ind for ind in toSolve if not name == transToSeq[ind][0]]

def solveFASTA(in_fasta, out_fasta, tid_regions, transToSeq):
    """
    Makes transcripts for every sequence in the input file.
    """
    header = ""
    seq = ""
    toSolve = list(tid_regions.keys())
    for line in in_fasta:
        line = line.strip()
        if len(line) == 0:
            continue
        if line.startswith(">"):
            if len(header) > 0:
                toSolve = solveSeq(out_fasta, header, seq, toSolve, tid_regions, transToSeq)
                seq = ""
            header = line
        else:
            seq = seq + line
    if len(header) > 0:
        toSolve = solveSeq(out_fasta, header, seq, toSolve, tid_regions, transToSeq)

def count(dict):
    n = 0
    for idx in dict:
        n = n + len(dict[idx])
    return n

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Error - 3 arguments needed: in.gtf in.fa out.fa [flag1, flag2...]"
        sys.exit(0)

    for arg in sys.argv[3::]:
        if arg == 's':
            LOOK_S = True
        elif arg == 'alt':
            ALT = True

    gff = open(sys.argv[1])
    tid_exons, transToSeq = parse(gff)
    gff.close()

    tid_regions = makeRegions(tid_exons)

    in_fasta = open(sys.argv[2])
    out_fasta = open(sys.argv[3],'w')
    solveFASTA(in_fasta, out_fasta, tid_regions, transToSeq)
    in_fasta.close()
    out_fasta.close()

    print "plus: " + str(plus)
    print "minus: " + str(minus)
    print "duljina: " + str(duljina)
    
