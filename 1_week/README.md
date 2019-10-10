# Algorithms for DNA Sequencing
## Week 1: DNA Sequencing, strings and matching
### Lectures
1. [Introduction](docs/intro.pdf)
2. [Why Study This?](docs/why_study_this.pdf)
3. [DNA sequencing past and present](docs/DNA_seq_past_present.pdf)
4. [Genomes as strings, reads as substrings](docs/genomes_strings.pdf)
5. [String definitions and Python examples](docs/python_str_def_ex.pdf)
6. [How DNA gets copied](docs/DNA_copying.pdf)
7. [How second-generation sequencers work](docs/second_gen_parallel.pdf)
8. [Sequencing errors and base qualities](docs/seq_errors_base_qualities.pdf)
9. [Working with sequencing reads (FASTQ format)](docs/fastq_format.pdf)
10. [Sequencers give pieces to genomic puzzles](docs/pieces_fragmentary.pdf)
    * Unrelated humans have genomes of ~3 billion base pairs which are 99.8 - 99.9% similar
11. [Read alignment and why it's hard](docs/read_alignment_hard.pdf)
12. [Naive exact matching](docs/naive_exact_matching.pdf)

### Resources
* [ERR037900_1.first1000.fastq](1_week/ERR037900_1.first1000.fastq)
* [lambda_virus.fa](1_week/lambda_virus.fa)

### Utility Functions
```python
def QtoPhred33(Q):
    """Turn Q into Phred+33 ASCII-­‐encoded quality"""
    return chr(Q + 33) # converts character to integer according to ASCII table

def phred33ToQ(qual):
    """Turn Phred+33 ASCII-encoded quality into Q"""
    return ord(qual) - 33 # converts integer to character according to ASCII table
```

### Assignment 1
In lecture and in a practical, we saw an implementation of the naive exact matching algorithm:
```python
def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences
```
...and we saw a function that takes a DNA string and returns its reverse complement:
```python
def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t
```

First, implement a version of the naive exact matching algorithm that is strand-aware. That is, instead of looking only for occurrences of P in T, additionally look for occurrences of the reverse complement of P in T. If P is ACT, your function should find occurrences of both ACTand its reverse complement AGT in T.

If P and its reverse complement are identical (e.g. AACGTT), then a given match offset should be reported only once. So if your new function is called naive_with_rc, then the old naive function and your new naive_with_rc function should return the same results when P equals its reverse complement.

Next, download and parse the [lambda virus genome](lambda_virus.fa).

### Solution 1
* [.ipynb](1_assignment.ipynb)
* [.py](1_assignment.py)
* [.html](https://nbviewer.jupyter.org/github/claytonjwong/Algorithms-DNA-Sequencing/blob/master/1_week/1_assignment.ipynb)