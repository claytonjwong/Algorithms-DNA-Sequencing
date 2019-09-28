#
# Practical: Implementating a K-mer Index
#
import bisect

class Index(object):
    def __init__(self, t, k):
        self.k = k
        self.index = []
        for i in range(len(t)-k+1):
            self.index.append((t[i:i+k], i))
        self.index.sort()

    def query(self, p):
        kmer = p[:self.k]
        i = bisect.bisect_left(self.index, (kmer, -1))
        hits = []
        while i < len(self.index):
            if self.index[i][0] != kmer:
                break # end of multimap equal range
            hits.append(self.index[i][1])
            i += 1
        return hits

def queryIndex(p, t, index):
    k = index.k
    offsets = []
    for i in index.query(p):
        if (p[k:] == t[i+k:i+len(p)]):
            offsets.append(i)
    return offsets

t = 'GCTACGATCTAGAATCTA'
p = 'TCTA'

index = Index(t, 2)
ans = queryIndex(p, t, index)
print(ans)

#
# partial matching algroithm implementation using an exact matching algorithm (Indexing)
# combined with the pigeon hole principle to allow up to k mismatches of pattern in text
#
def approximate_match_index(p, t, k):
    segment_length = round(len(p) // (k+1))
    all_matches = set()
    index = Index(t, 8) # built on 8-mers
    for i in range(k+1):
        start = i * segment_length
        end = min((i+1) * segment_length, len(p))
        matches = index.query(p[start:end])
        for m in matches:
            text_offset = m - start
            if text_offset < 0 or (text_offset + len(p)) > len(t):
                continue
            mismatches = 0
            for j in range(0, start):
                if not p[j] == t[text_offset + j]:
                    mismatches += 1
                    if mismatches > k:
                        break
            for j in range(end, len(p)):
                if not p[j] == t[text_offset + j]:
                    mismatches += 1
                    if mismatches > k:
                        break
            if mismatches <= k:
                all_matches.add(text_offset)
    return list(all_matches)

def read_FAST_A(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


p = "GGCGCGGTGGCTCACGCCTGTAAT"    
t = read_FAST_A('chr1.GRCh38.excerpt.fasta')
print(approximate_match_index(p, t, 2))
    


class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits
    