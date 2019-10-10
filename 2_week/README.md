# Algorithms for DNA Sequencing
## Week 2: Preprocessing, Indexing, and Approximate Matching
### Lectures
1. [Boyer-Moore: Basics](docs/boyer-moore.pdf)
2. [Boyer-Moore: Putting It All Together](docs/boyer-moore-together.pdf)
3. [Diversion: Repetitive Elements](docs/repetitive_elements.pdf)
4. [Preprocessing](docs/preprocessing.pdf)
5. [Indexing and K-mers](docs/indexing_kmers.pdf)
6. [Data Structures for Indexing](docs/data_structures.pdf)
7. [Hash Tables for Indexing](docs/hash_tables.pdf)
8. [Variations on K-mer Indexes](docs/indexing_variations.pdf)
9. [Indexing by Suffix](docs/suffix.pdf)
10. [Approximate Matching, Hamming, and Edit Distance](docs/approximate.pdf)
11. [Pigeonhole Principle](docs/pigeonhole.pdf)

### Resources
* [bm_preproc.py](bm_preproc.py)
* [k-mer_index.py](k-mer_index.py)
* [chr1.GRCh38.excerpt.fasta](chr1.GRCh38.excerpt.fasta)

### Assignment 2
In a practical, we saw Python code implementing the Boyer-Moore algorithm. Some of the code is for preprocessing the pattern P into the tables needed to execute the bad character and good suffix rules — we did not discuss that code. But we did discuss the code that performs the algorithm given those tables:
```python
def boyer_moore(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
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
    return occurrences
```

**Measuring Boyer-Moore's benefit.** First, download the Python module for [Boyer-Moore preprocessing](bm_preproc.py)

This module provides the BoyerMoore class, which encapsulates the preprocessing info used by the boyer_moore function above. Second, download the provided excerpt of [human chromosome 1](chr1.GRCh38.excerpt.fasta).

Third, implement versions of the naive exact matching and Boyer-Moore algorithms *that additionally count and return (a) the number of character comparisons performed and (b) the number of alignments tried.* Roughly speaking, these measure how much work the two different algorithms are doing.

#### Naive
```python
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
```

#### Boyer-Moore
```python
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
```

#### Approximate Matching with Boyer Moore
```python
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
```

#### Approximate Matching with Indexing
```python
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
```

#### Approximate Matching with Subsequence Indexing
```python
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
```

### Solution 2
* [.ipynb](2_assignment.ipynb)
* [.py](2_assignment.py)
* [.html](https://nbviewer.jupyter.org/github/claytonjwong/Algorithms-DNA-Sequencing/blob/master/2_week/2_assignment.ipynb)
