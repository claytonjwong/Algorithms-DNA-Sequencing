#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python

"""bm_preproc.py: Boyer-Moore preprocessing."""

__author__ = "Ben Langmead"

import unittest


def z_array(s):
    """ Use Z algorithm (Gusfield theorem 1.4.1) to preprocess s """
    assert len(s) > 1
    z = [len(s)] + [0] * (len(s)-1)

    # Initial comparison of s[1:] with prefix
    for i in range(1, len(s)):
        if s[i] == s[i-1]:
            z[1] += 1
        else:
            break

    r, l = 0, 0
    if z[1] > 0:
        r, l = z[1], 1

    for k in range(2, len(s)):
        assert z[k] == 0
        if k > r:
            # Case 1
            for i in range(k, len(s)):
                if s[i] == s[i-k]:
                    z[k] += 1
                else:
                    break
            r, l = k + z[k] - 1, k
        else:
            # Case 2
            # Calculate length of beta
            nbeta = r - k + 1
            zkp = z[k - l]
            if nbeta > zkp:
                # Case 2a: zkp wins
                z[k] = zkp
            else:
                # Case 2b: Compare characters just past r
                nmatch = 0
                for i in range(r+1, len(s)):
                    if s[i] == s[i - k]:
                        nmatch += 1
                    else:
                        break
                l, r = k, r + nmatch
                z[k] = r - k + 1
    return z


def n_array(s):
    """ Compile the N array (Gusfield theorem 2.2.2) from the Z array """
    return z_array(s[::-1])[::-1]


def big_l_prime_array(p, n):
    """ Compile L' array (Gusfield theorem 2.2.2) using p and N array.
        L'[i] = largest index j less than n such that N[j] = |P[i:]| """
    lp = [0] * len(p)
    for j in range(len(p)-1):
        i = len(p) - n[j]
        if i < len(p):
            lp[i] = j + 1
    return lp


def big_l_array(p, lp):
    """ Compile L array (Gusfield theorem 2.2.2) using p and L' array.
        L[i] = largest index j less than n such that N[j] >= |P[i:]| """
    l = [0] * len(p)
    l[1] = lp[1]
    for i in range(2, len(p)):
        l[i] = max(l[i-1], lp[i])
    return l


def small_l_prime_array(n):
    """ Compile lp' array (Gusfield theorem 2.2.4) using N array. """
    small_lp = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i+1:  # prefix matching a suffix
            small_lp[len(n)-i-1] = i+1
    for i in range(len(n)-2, -1, -1):  # "smear" them out to the left
        if small_lp[i] == 0:
            small_lp[i] = small_lp[i+1]
    return small_lp


def good_suffix_table(p):
    """ Return tables needed to apply good suffix rule. """
    n = n_array(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)


def good_suffix_mismatch(i, big_l_prime, small_l_prime):
    """ Given a mismatch at offset i, and given L/L' and l' arrays,
        return amount to shift as determined by good suffix rule. """
    length = len(big_l_prime)
    assert i < length
    if i == length - 1:
        return 0
    i += 1  # i points to leftmost matching position of P
    if big_l_prime[i] > 0:
        return length - big_l_prime[i]
    return length - small_l_prime[i]


def good_suffix_match(small_l_prime):
    """ Given a full match of P to T, return amount to shift as
        determined by good suffix rule. """
    return len(small_l_prime) - small_l_prime[1]


def dense_bad_char_tab(p, amap):
    """ Given pattern string and list with ordered alphabet characters, create
        and return a dense bad character table.  Table is indexed by offset
        then by character. """
    tab = []
    nxt = [0] * len(amap)
    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i+1
    return tab


class BoyerMoore(object):
    """ Encapsulates pattern and associated Boyer-Moore preprocessing. """

    def __init__(self, p, alphabet='ACGT'):
        # Create map from alphabet characters to integers
        self.amap = {alphabet[i]: i for i in range(len(alphabet))}
        # Make bad character rule table
        self.bad_char = dense_bad_char_tab(p, self.amap)
        # Create good suffix rule table
        _, self.big_l, self.small_l_prime = good_suffix_table(p)

    def bad_character_rule(self, i, c):
        """ Return # skips given by bad character rule at offset i """
        assert c in self.amap
        assert i < len(self.bad_char)
        ci = self.amap[c]
        return i - (self.bad_char[i][ci]-1)

    def good_suffix_rule(self, i):
        """ Given a mismatch at offset i, return amount to shift
            as determined by (weak) good suffix rule. """
        length = len(self.big_l)
        assert i < length
        if i == length - 1:
            return 0
        i += 1  # i points to leftmost matching position of P
        if self.big_l[i] > 0:
            return length - self.big_l[i]
        return length - self.small_l_prime[i]

    def match_skip(self):
        """ Return amount to shift in case where P matches T """
        return len(self.small_l_prime) - self.small_l_prime[1]


class TestBoyerMoorePreproc(unittest.TestCase):

    def test_z_1(self):
        s = 'abb'
        #    -00
        z = z_array(s)
        self.assertEqual([3, 0, 0], z)

    def test_z_2(self):
        s = 'abababab'
        #    00604020
        z = z_array(s)
        self.assertEqual([8, 0, 6, 0, 4, 0, 2, 0], z)

    def test_z_3(self):
        s = 'abababab'
        #    00604020
        z = z_array(s)
        self.assertEqual([8, 0, 6, 0, 4, 0, 2, 0], z)

    def test_n_1(self):
        s = 'abb'
        #    01-
        n = n_array(s)
        self.assertEqual([0, 1, 3], n)

    def test_n_2(self):
        s = 'abracadabra'
        #    1004010100-
        n = n_array(s)
        self.assertEqual([1, 0, 0, 4, 0, 1, 0, 1, 0, 0, 11], n)

    def test_n_3(self):
        s = 'abababab'
        #    0204060-
        n = n_array(s)
        self.assertEqual([0, 2, 0, 4, 0, 6, 0, 8], n)

    def test_big_l_prime_1(self):
        s = 'abb'
        #    001
        big_l_prime = big_l_prime_array(s, n_array(s))
        self.assertEqual([0, 0, 2], big_l_prime)

    def test_big_l_prime_2(self):
        s = 'abracadabra'
        #    01234567890
        # L' 00000003007
        # L  00000003337
        big_l_prime = big_l_prime_array(s, n_array(s))
        self.assertEqual([0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 8], big_l_prime)

    def test_small_l_prime_1(self):
        s = 'abracadabra'
        # N  1004010100-
        # l'           1
        # l'        4
        # l' 44444444111
        small_l_prime = small_l_prime_array(n_array(s))
        self.assertEqual([11, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1], small_l_prime)

    def test_good_suffix_match_mismatch_1(self):
        p = 'GGTAGGT'
        big_l_prime, big_l, small_l_prime = good_suffix_table(p)
        self.assertEqual([0, 0, 0, 0, 3, 0, 0], big_l_prime)
        self.assertEqual([0, 0, 0, 0, 3, 3, 3], big_l)
        self.assertEqual([7, 3, 3, 3, 3, 0, 0], small_l_prime)
        self.assertEqual(0, good_suffix_mismatch(6, big_l_prime, small_l_prime))
        self.assertEqual(0, good_suffix_mismatch(6, big_l, small_l_prime))
        #  t:      xT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(7, good_suffix_mismatch(5, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(5, big_l, small_l_prime))
        #  t:     xGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(7, good_suffix_mismatch(4, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(4, big_l, small_l_prime))
        #  t:    xGGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(4, good_suffix_mismatch(3, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(3, big_l, small_l_prime))
        #  t:   xAGGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(4, good_suffix_mismatch(2, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(2, big_l, small_l_prime))
        #  t:  xTAGGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(4, good_suffix_mismatch(1, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(1, big_l, small_l_prime))
        #  t: xGTAGGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(4, good_suffix_mismatch(0, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(0, big_l, small_l_prime))

    def test_good_suffix_table_1(self):
        s = 'abb'
        #    001
        big_l_prime, big_l, small_l_prime = good_suffix_table(s)
        self.assertEqual([0, 0, 2], big_l_prime)
        self.assertEqual([0, 0, 2], big_l)
        self.assertEqual([3, 0, 0], small_l_prime)

    def test_good_suffix_table_2(self):
        s = 'abracadabra'
        #    01234567890
        # L' 00000003007
        # L  00000003337
        # l' -4444444111
        big_l_prime, big_l, small_l_prime = good_suffix_table(s)
        self.assertEqual([0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 8], big_l_prime)
        self.assertEqual([0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 8], big_l)
        self.assertEqual([11, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1], small_l_prime)
        


# In[2]:


def boyer_moore(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    alignments = 0
    comparisons = 0
    while i < len(t) - len(p) + 1:
        alignments += 1
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            comparisons += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences, alignments, comparisons


# In[3]:


def naive(p, t):
    occurrences = []
    alignments = 0
    comparisons = 0
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        alignments += 1
        match = True
        for j in range(len(p)):  # loop over characters
            comparisons += 1
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences, alignments, comparisons


# In[4]:


def read_FAST_A(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


# In[5]:


human_chromosome_1 = read_FAST_A('chr1.GRCh38.excerpt.fasta')


# In[6]:


print(human_chromosome_1)


# In[7]:


occurrences, alignments, comparisons = naive('GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG', human_chromosome_1)


# In[8]:


question1 = '''How many alignments does the naive exact matching algorithm try when matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to the excerpt of human chromosome 1? (Don't consider reverse complements.)'''


# In[9]:


print(question1)


# In[10]:


print(alignments)


# In[11]:


question2 = '''How many character comparisons does the naive exact matching algorithm try when matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to the excerpt of human chromosome 1? (Don't consider reverse complements.)'''


# In[12]:


print(question2)


# In[13]:


print(comparisons)


# In[14]:


question3 = '''How many alignments does Boyer-Moore try when matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to the excerpt of human chromosome 1? (Don't consider reverse complements.)'''


# In[15]:


print(question3)


# In[16]:


occurrences, alignments, comparisons = boyer_moore('GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG', BoyerMoore('GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'), human_chromosome_1)


# In[17]:


print(alignments)


# In[74]:


#
# Practical: partial matching algroithm implementation using an exact matching algorithm (Boyer-Moore)
# combined with the pigeon hole principle to allow up to k mismatches of pattern in text
#
def approximate_match_boyer_moore(p, t, k):
    segment_length = round(len(p) // (k+1))
    all_matches = set()
    for i in range(k+1):
        start = i * segment_length
        end = min((i+1) * segment_length, len(p))
        p_bm = BoyerMoore(p[start:end], alphabet='ACGT')
        matches, _, _ = boyer_moore(p[start:end], p_bm, t)
        for m in matches:
            #
            # Question: why is text offset equal to the matched position minus start of current segment?
            #
            # Answer: this is beacuse 'start' is the k-th segment offset within p, so we are basically aligning pattern within text based on the k-th exact matching segment
            # so text offset from T's perspective is where P would begin to align within T.
            #
            # Example:
            #
            #       01234567
            #   T = ABCDEFG
            #   P =   CDEF
            #   K = 1
            #
            # since K = 1, then K + 1 = 2, thus P is split into 2 partitions CD and EF
            # the first partition:  CD matches T at m = 2, the start of partition CD is 0, so P aligns with T at m = 2 - 0 = 2
            # the second partition: EF matches T at m = 4, the start of partition EF is 2, so P aligns with T at m = 4 - 2 = 2
            # thus for each match in a k-th partition, P aligns at same offset m within T ( ie. offset 2 for this example )
            #
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


# In[56]:


p = 'AACTTG'
t = 'CACTTAATTTG'
print(approximate_match_boyer_moore(p, t, 2))


# In[57]:


#
# Practical: Implementating a K-mer Index
#
"""kmer_index.py: A k-mer index for indexing a text."""

__author__ = "Ben Langmead"

import bisect

class Index(object):
    """ Holds a substring index for a text T """
    def __init__(self, t, k):
        """ Create index from all substrings of t of length k """
        self.k = k # k-mer length (k)
        self.index = []
        for i in range(len(t)-k+1): # for each k-mer
            self.index.append((t[i:i+k], i)) # add (k-mer, offset) pair
        self.index.sort()

    def query(self, p):
        """ Return index hits for first k-mer of p """
        kmer = p[:self.k] # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1)) # binary search
        hits = []
        while i < len(self.index): # # collect matching index entries
            if self.index[i][0] != kmer:
                break # end of multimap equal range
            hits.append(self.index[i][1])
            i += 1
        return hits


# In[58]:


def queryIndex(p, t, index):
    k = index.k
    offsets = []
    for i in index.query(p):
        if (p[k:] == t[i+k:i+len(p)]):
            offsets.append(i)
    return offsets


# In[59]:


t = 'GCTACGATCTAGAATCTA'
p = 'TCTA'


# In[60]:


index = Index(t, 2)
ans = queryIndex(p, t, index)
print(ans)


# In[61]:


print(t[7:7+len(p)])


# In[62]:


print(t[14:14+len(p)])


# In[72]:


question4_part1 = '''Implement the pigeonhole principle using Index to find exact matches for the partitions.
Assume P always has length 24, and that we are looking for approximate matches with up to 2 mismatches
(substitutions). We will use an 8-mer index.  Write a function that, given a length-24 pattern P
and given an Index object built on 8-mers, finds all approximate occurrences of P within T with up to 2 mismatches.
Insertions and deletions are not allowed. Don't consider any reverse complements.'''


# In[73]:


print(question4_part1)


# In[101]:


#
# partial matching algroithm implementation using an exact matching algorithm (Indexing)
# combined with the pigeon hole principle to allow up to k mismatches of pattern in text
#
def approximate_match_index(p, t, k):
    segment_length = round(len(p) // (k+1))
    hits = 0
    all_matches = set()
    index = Index(t, 8) # built on 8-mers
    for i in range(k+1):
        start = i * segment_length
        end = min((i+1) * segment_length, len(p))
        matches = index.query(p[start:end])
        hits += len(matches)
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
    return list(all_matches), hits


# In[102]:


question4_part2 = '''How many times does the string GGCGCGGTGGCTCACGCCTGTAAT, which is derived
from a human Alu sequence, occur with up to 2 substitutions in the excerpt of human chromosome 1?
(Don't consider reverse complements here.)'''


# In[103]:


print(question4_part2)


# In[109]:


t = human_chromosome_1
p = 'GGCGCGGTGGCTCACGCCTGTAAT'
k = 2
ans, hits = approximate_match_index(p, t, k)
print("occurrences: " + str(len(ans)))


# In[106]:


question5 = '''Using the instructions given in Question 4, how many total index hits are there
when searching for occurrences of GGCGCGGTGGCTCACGCCTGTAAT with up to 2 substitutions in the excerpt
of human chromosome 1?

(Don't consider reverse complements.)

Hint: You should be able to use the boyer_moore function (or the slower naive function)
to double-check your answer.'''


# In[107]:


print(question5)


# In[110]:


print("index hits: " + str(hits))


# In[112]:


question6_part1 = '''Let's examine whether there is a benefit to using an index built using subsequences of T
rather than substrings, as we discussed in the "Variations on k-mer indexes" video. We'll consider subsequences
involving every N characters. For example, if we split ATATAT into two substring partitions, we would get partitions
ATA (the first half) and TAT (second half). But if we split ATATAT into two subsequences by taking every other
character, we would get AAA (first, third and fifth characters) and TTT (second, fourth and sixth).

Another way to visualize this is using numbers to show how each character of P is allocated to a partition.
Splitting a length-6 pattern into two substrings could be represented as 111222, and splitting into two subsequences
of every other character could be represented as 121212

The following class SubseqIndex is a more general implementation of Index that additionally handles subsequences.
It only considers subsequences that take every Nth character:
'''


# In[113]:


print(question6_part1)


# In[114]:


import bisect
   
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


# In[115]:


ind = SubseqIndex('ATATAT', 3, 2)
print(ind.index)


# In[116]:


p = 'TTATAT' 
print(ind.query(p[0:])) # TAA (ie. each even index character is NOT in the index)


# In[118]:


p = 'TTATAT' 
print(ind.query(p[1:])) # TTT (ie. each odd index character is in the index)


# In[128]:


question6_part2 = '''Write a function that, given a length-24 pattern P and given a SubseqIndex object
built with k = 8 and ival = 3, finds all approximate occurrences of P within T with up to 2 mismatches.

When using this function, how many total index hits are there when searching for GGCGCGGTGGCTCACGCCTGTAAT
with up to 2 substitutions in the excerpt of human chromosome 1? (Again, don't consider reverse complements.)'''


# In[129]:


print(question6_part2)


# In[126]:


#
# partial matching algroithm implementation using an exact matching algorithm (Subsequence Indexing)
# combined with the pigeon hole principle to allow up to k mismatches of pattern in text
#
# note: in this solution, the partitions overlap ( ie. every 3rd subsequence interval, 123123123 )
# compared to the previous solutions where each partition has a unique, non-overlapping start/end
# ( ie. every 3rd chunk, 111222333 )
#
def approximate_match_subseq_index(p, t, k):
    segment_length = round(len(p) // (k+1))
    hits = 0
    all_matches = set()
    index = SubseqIndex(t, 8, 3) # built on 8-mers and subsequence intervals of 3
    for start in range(k+1):
        matches = index.query(p[start:])
        hits += len(matches)
        for m in matches:
            text_offset = m - start
            if text_offset < 0 or (text_offset + len(p)) > len(t):
                continue
            mismatches = 0
            for j in range(0, len(p)):
                if not p[j] == t[text_offset + j]:
                    mismatches += 1
                    if mismatches > k:
                        break
            if mismatches <= k:
                all_matches.add(text_offset)
    return list(all_matches), hits


# In[127]:


t = human_chromosome_1
p = 'GGCGCGGTGGCTCACGCCTGTAAT'
k = 2
_, hits = approximate_match_subseq_index(p, t, k)
print(hits)


# In[ ]:




