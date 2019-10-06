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

### Solution 2
* [.ipynb](2_assignment.ipynb)
* [.py](2_assignment.py)
* [.md](2_assignment.md)
