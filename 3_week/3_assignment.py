#!/usr/bin/env python
# coding: utf-8

# In[8]:


print('''I wrote this code de novo (ie. from scratch).  I feel this code is easier to read and understand
than the lecture's implementation.  I also wrote a memo for the top-down solution to run efficiently
without repeatedly solving overlapping sub-problems.

Let i, j be the index one-past the end of the prefix substring (ie. 0..i-1 inclusive and 0..j-1 inclusive)

Strings of length m and n can be denoted as 0..i and 0..j initially when i = m and j = n

*example:*

    A = abc
        0123
           m
           i

    A[0..i) == A[0..i-1] == abc 
               
    B = wxyz
        01234
            n
            j
    
    B[0..j) == B[0..j-1] == wxyz
    
Assume we know the optimal solution for the prefix of the strings A, B without considering the last character
of each string.  Then take the last character of each string into consideration.  The last character
from each string can be either a substitution (match/mismatch), insertion, or deletion.

The base case occurs when either A or B is an empty string.  The distance from an empty string p
to another non-empty string q is the length of q.''')


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


# In[9]:


print('''
### Assignment 3
#### Questions 1 and 2
We saw how to adapt dynamic programming to find approximate occurrences of a pattern in a text. Recall that:
* Rows of the dynamic programming matrix are labeled with bases from P and columns with bases from T
* Elements in the **first row are set to 0**
* Elements in the first column are set to 0, 1, 2, ..., (same as edit distance)
* Other elements are set in the same way as elements of a standard edit distance matrix
* The **minimal value in the bottom row** is the edit distance of the closest match between P and T

First, download the provided excerpt of [human chromosome 1](chr1.GRCh38.excerpt.fasta)

Second, parse it using the readGenome function we wrote before.

Third, adapt the editDistance function we saw in practical to answer questions 1 and 2 below.
Your function should take arguments p (pattern), t (text) and should return the edit distance
of the match between P and T with the fewest edits.''')


# In[63]:


get_ipython().run_cell_magic('time', '', 'print(editDistanceNaive("cha-ching", "cha-ching!!!"))')


# In[64]:


get_ipython().run_cell_magic('time', '', 'print(editDistanceMemo("cha-ching", "cha-ching!!!"))')


# In[65]:


get_ipython().run_cell_magic('time', '', 'print(editDistanceDP("cha-ching", "cha-ching!!!"))')


# In[4]:


print('''Primary difference between edit distance and global alignment is the ability
to apply a weighted penalty to substitutions and indels (insertion/deletions)
''')


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


# In[10]:


print('''
Say we are concerned only with overlaps that (a) are exact matches (no differences allowed),
and (b) are at least K bases long. To make an overlap graph, we could call 
overlap(a,b,min_length=k) on every possible pair of reads from the dataset.
Unfortunately, that will be very slow!

Consider this: Say we are using k=6, and we have a read A whose length-6 suffix is GTCCTA.
Say GTCCTA does not occur in any other read in the dataset. In other words, the 6-mer GTCCTA
occurs at the end of read A and nowhere else. It follows that A's suffix cannot possibly overlap
the prefix of any other read by 6 or more characters.

Put another way, if we want to find the overlaps involving a suffix of read A and a prefix of some other read,
we can ignore any reads that don't contain the length-k suffix of A. This is good news because it can save us
a lot of work!

Here is a suggestion for how to implement this idea. You don't have to do it this way, but this might help you.
Let every k-mer in the dataset have an associated Python Set object, which starts out empty.  We use a Python dictionary
to associate each k-mer with its corresponding Set.

1. For every k-mer in a read, we add the read
to the Set object corresponding to that k-mer. If our read is GATTA and k=3, we would add GATTA to the Set objects
for GAT, ATT and TTA. We do this for every read so that, at the end, each Set contains all reads containing
the corresponding k-mer.
2. Now, for each read A, we find all overlaps involving a suffix of A. To do this, we take A's length-k suffix,
find all reads containing that k-mer (obtained from the corresponding Set) and call overlap(a,b,min_length=k)
for each.

The most important point is that we do not call overlap(a,b,min_length=k) if B does not contain
the length-k suffix of A.

Download and parse the read sequences from the provided Phi-X FASTQ file. We'll just use their base sequences,
so you can ignore read names and base qualities. Also, no two reads in the FASTQ have the same sequence of bases.
This makes things simpler.

Next, find all pairs of reads with an exact suffix/prefix match of length at least 30. Don't overlap a read
with itself; if a read has a suffix/prefix match to itself, ignore that match. Ignore reverse complements.
* Hint 1: Your function should not take much more than 15 seconds to run on this 10,000-read dataset,
and maybe much less than that. (Our solution takes about 3 seconds.) If your function is much slower,
there is a problem somewhere.
* Hint 2: Remember not to overlap a read with itself. If you do, your answers will be too high.
* Hint 3: You can test your implementation by making up small examples, then checking that
(a) your implementation runs quickly, and (b) you get the same answer as if you had simply called
overlap(a,b,min_length=k) on every pair of reads.  We also have provided a couple
[examples you can check against](https://nbviewer.jupyter.org/github/BenLangmead/ads1-hw-examples/blob/master/hw3_overlap_all.ipynb).''')


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


# In[11]:


print('''
* **Question 3:** Picture the overlap graph corresponding to the overlaps just calculated. How many edges are in the graph?
In other words, how many distinct pairs of reads overlap?

* **Question 4:** Picture the overlap graph corresponding to the overlaps computed for the previous question.
How many nodes in this graph have at least one outgoing edge? (In other words, how many reads have a suffix
involved in an overlap?)''')


# In[51]:


reads, _ = readFAST_Q('ERR266411_1.for_asm.fastq')
pairs = overlap_all_pairs(reads, 30)


# In[52]:


print(len(pairs))


# In[53]:


print(len(set(pair[0] for pair in pairs)))


# In[ ]:




