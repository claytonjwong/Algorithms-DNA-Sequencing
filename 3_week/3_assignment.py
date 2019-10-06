#!/usr/bin/env python
# coding: utf-8

# In[12]:


def editDistanceNaive(A, B):
    def go(A, B, i, j):
        if i == 0: return j
        if j == 0: return i
        return min(
            go(A, B, i-1, j-1) + int(A[i-1] != B[j-1]),
            go(A, B, i, j-1) + 1,
            go(A, B, i-1, j) + 1,
        )
    return go(A, B, len(A), len(B))

def editDistanceMemo(A, B):
    def go(A, B, i, j, memo={}):
        key = str(i) + ',' + str(j)
        if memo.get(key): return memo[key]
        elif i == 0: memo[key] = j
        elif j == 0: memo[key] = i
        else:
            memo[key] = min(
                go(A, B, i-1, j-1, memo) + int(A[i-1] != B[j-1]),
                go(A, B, i, j-1, memo) + 1,
                go(A, B, i-1, j, memo) + 1,
            )
        return memo[key]
    return go(A, B, len(A), len(B))

def editDistanceDP(A, B):
    m, n = len(A), len(B)
    dp = [[0 for j in range(n+1)] for i in range(m+1)]
    for i in range(m+1): dp[i][0] = i
    for j in range(n+1): dp[0][j] = j
    for i in range(1, m+1):
        for j in range(1, n+1):
            dp[i][j] = min(
                dp[i-1][j-1] + int(A[i-1] != B[j-1]),
                dp[i-1][j] + 1,
                dp[i][j-1] + 1,
            )
    return dp[m][n]


# In[63]:


get_ipython().run_cell_magic('time', '', 'print(editDistanceNaive("cha-ching", "cha-ching!!!"))')


# In[64]:


get_ipython().run_cell_magic('time', '', 'print(editDistanceMemo("cha-ching", "cha-ching!!!"))')


# In[65]:


get_ipython().run_cell_magic('time', '', 'print(editDistanceDP("cha-ching", "cha-ching!!!"))')


# In[69]:


def editDistanceApproximate(P, T):
    m, n = len(P), len(T)
    dp = [[0 for j in range(n+1)] for i in range(m+1)]
    for i in range(m+1): dp[i][0] = i # init first column by distance from empty string
        
#   the first row is all 0s unlike edit distance, since there is no bias toward alignment
#   of P in T from the beginning of both P and T, (ie. P can start at any index in T)

#
#   for j in range(n+1): dp[0][j] = j # DELETED!!!
#

    for i in range(1, m+1):
        for j in range(1, n+1):
            dp[i][j] = min(
                dp[i-1][j-1] + int(P[i-1] != T[j-1]),
                dp[i-1][j] + 1,
                dp[i][j-1] + 1,
            )
    return min(dp[m])


# In[70]:


P = "GCGTATGC"
T = "TATTGGCTATACGGTT"
print(editDistanceApproximate(P, T))


# In[71]:


print('First, download the provided excerpt of human chromosome 1.  Second, parse it using the readGenome function we wrote before.')


# In[72]:


def read_FAST_A(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


# In[73]:


T = read_FAST_A("chr1.GRCh38.excerpt.fasta")


# In[74]:


print(T)


# In[75]:


question1 = 'What is the edit distance of the best match between pattern GCTGATCGATCGTACG and the excerpt of human chromosome 1? (Don\'t consider reverse complements.)'
print(question1)


# In[76]:


print(editDistanceDP("GCTGATCGATCGTACG", T))


# In[77]:


question2 = 'What is the edit distance of the best match between pattern GATTTACCAGATTGAG and the excerpt of human chromosome 1? (Don\'t consider reverse complements.)'


# In[78]:


print(editDistanceDP("GATTTACCAGATTGAG", T))


# In[113]:


question3_part1 = '''In a practical, we saw a function for finding the longest exact overlap (suffix/prefix match)
between two strings. The function is copied below.'''
print(question3_part1)


# In[29]:


def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match


# In[30]:


def readFAST_Q(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


# In[43]:


reads, _ = readFAST_Q('ERR266411_1.for_asm.fastq')


# In[44]:


print(reads)


# In[45]:


def overlap_all_pairs(reads, k, map={}):
    def get_kmers(read, k):
        res = set()
        for i in range(0, len(read)-k+1):
            res.add(read[i:i+k])
        return res
    for read in reads:
        kmers = get_kmers(read, k)
        for kmer in kmers:
            if not kmer in map.keys():
                map[kmer] = set()
            map[kmer].add(read)
    pairs = []
    for head in reads:
        kmer = head[-k:]
        candidates = map[kmer]
        for tail in candidates:
            if (not head == tail and overlap(head, tail, k)):
                pairs.append((head, tail))
    return pairs


# In[46]:


print(overlap_all_pairs(['ABCDEFG', 'EFGHIJ', 'HIJABC'], 3))


# In[47]:


print(overlap_all_pairs(['ABCDEFG', 'EFGHIJ', 'HIJABC'], 4))


# In[50]:


print(overlap_all_pairs(['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT'], 4))


# In[51]:


reads, _ = readFAST_Q('ERR266411_1.for_asm.fastq')
pairs = overlap_all_pairs(reads, 30)


# In[52]:


print(len(pairs))


# In[53]:


print(len(set(pair[0] for pair in pairs)))


# In[ ]:




