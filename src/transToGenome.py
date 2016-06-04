import utility_sam
import copy
import sys
import parser
import re;

def checkBitMask(name, bitMask, tid_regions):
    bmLen = sum([1 for c in bitMask if c == '1' or c == '0'])
    if bmLen == len(bitMask) and bmLen == len(tid_regions[name]):
        return True
    return False

def extractRegions2(tid, tid_regions):
    fields = tid.split("_")
    if len(fields) != 1:
        name = "_".join(fields[0:len(fields)-1])
        bm = fields[-1]
        if checkBitMask(name, bm, tid_regions):
            return name, [region for i, region in enumerate(tid_regions[name]) if bm[i] == '1']
    return tid, tid_regions[tid]

def extractRegions(tid, tid_regions):
    if "_" not in tid:
        return tid, tid_regions[tid]
    name,bm = tid.split("_")
    return name, [region for i, region in enumerate(tid_regions[name]) if bm[i] == '1']

def solve(sam_lines, transToSeq, tid_regions):
    newSam = []
    for line in sam_lines:
        if line.rname == "*":
            newSam.append(line)
            continue

        name, regions = extractRegions(line.rname, tid_regions)
        cigar = line.SplitCigar()
        cigar.reverse() #stack
        regInd = 0
        curr = line.clipped_pos
        #In which region does it start
        last = 0
        while last + regions[regInd][1] - regions[regInd][0] + 1 < curr:
            last += regions[regInd][1] - regions[regInd][0] + 1
            regInd = regInd + 1

        posOnRef = regions[regInd][0] + curr - 1

        newCigar = ""
        regionSize = regions[regInd][1] - regions[regInd][0] + 1

        # last = 0
        while cigar:
            c,op = cigar.pop()
            count = int(c)
            if count + curr > regionSize + last:
                take = min(count, regionSize + last - curr + 1)
                #stavi intron
                if regInd < len(regions) - 1:
                    if take > 0:
                        newCigar += str(take) + op
                    if count - take > 0:
                        cigar.append((count - take, op))
                    newCigar += str(regions[regInd+1][0] - regions[regInd][1] - 1) + 'N'
                else:
                    newCigar += str(c) + op
                    break

                regInd = regInd + 1
                last += regionSize
                regionSize = regions[regInd][1] - regions[regInd][0] + 1
                if op != 'I':
                    curr += take
            else:
                if op != 'I':
                    curr += count
                newCigar += str(count) + op

        while cigar:
            c,op = cigar.pop()
            newCigar += str(c) + op

        newLine = copy.deepcopy(line);
        newLine.cigar = newCigar
        newLine.rname = transToSeq[name][0]
        m_front = re.match("^([\d]+)([SH])", newCigar)
        newLine.pos = posOnRef
        if m_front:
            newLine.pos += int(m_front.group(1))
        if transToSeq[name][1] == '-':
            newLine.flag ^= 0x10
        newSam.append(newLine)

    return newSam

if __name__ == "__main__":

    if len(sys.argv) < 4:
        print "Error - 3 arguments needed: in.sam in.gtf out.sam"
        sys.exit(0)

    headers, sam_lines = utility_sam.LoadSAM(sys.argv[1])
    gtf = open(sys.argv[2])
    tid_exons, transToSeq = parser.parse(gtf)
    gtf.close()
    tid_regions = parser.makeRegions(tid_exons)

    newSam = solve(sam_lines, transToSeq, tid_regions)
    samOut = open(sys.argv[3], "w")
    for head in headers:
        samOut.write(head + '\n')
    for sam in newSam:
        samOut.write(" ".join([sam.qname, str(sam.flag), sam.rname, str(sam.pos), \
            str(sam.mapq), sam.cigar, sam.mrnm, str(sam.mpos), str(sam.isize), sam.seq, sam.qual]) + "\n");
    samOut.close();
