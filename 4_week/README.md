# Algorithms for DNA Sequencing
## Week 4: Algorithms for Assembly
### Lectures
1. [Shortest Common Substring](4_week/docs/scss.pdf)
2. [Greedy Shortest Common Substring](4_week/docs/greedy_scss.pdf)
3. [3rd Law of Assembly: Repeats are Bad](4_week/docs/third_law.pdf)
4. [De Bruijn Graphs and Eulerian Walks](4_week/docs/dbg1.pdf)
5. [When Eulerian Walks Go Wrong](4_week/docs/dbg2.pdf)
6. [Assemblers in Practice](4_week/docs/assemblers_in_practice.pdf)
7. [The Future is Long?](4_week/docs/longreads.pdf)
8. [Wrap Up](4_week/docs/wrap_up.pdf)
    * Computer Scientists are signficiant contributors to Life Sciences since DNA, RNA, proteins, etc
    can be abstracted into algorithms based on strings.
    
### Assignment 4
#### Questions 1 and 2
In a practical, we saw the scs function (copied below along with overlap) for finding
the shortest common superstring of a set of strings.

```python
def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

import itertools

def scss(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest    
```

It's possible for there to be multiple different shortest common superstrings for the same set of input strings.
Consider the input strings ABC, BCA, CAB. One shortest common superstring is ABCAB but another is BCABC
and another is CABCA.

* **Question 1:** What is the length of the shortest common superstring of the following strings?
    * CCT, CTT, TGC, TGG, GAT, ATT

* **Question 2:** How many different shortest common superstrings are there for the input strings
given in the previous question?
    * Hint 1: You can modify the scs function to keep track of this.
    * Hint 2: You can look at [these examples](hw4_scss.ipynb) to double-check that your modified scs
    is working as expected.
    
#### Questions 3 and 4
Download this FASTQ file containing synthetic sequencing reads from a [mystery virus](mystery.fq).
All the reads are the same length (100 bases) and are exact copies of substrings from the forward strand
of the virus genome. You don't have to worry about sequencing errors, ploidy, or reads coming
from the reverse strand.

Assemble these reads using one of the approaches discussed, such as greedy shortest common superstring.
Since there are many reads, you might consider ways to make the algorithm faster, such as the one discussed
in the programming assignment in the previous module.

```python
def pick_max_overlap(reads, k):
    best_a, best_b, best_len = None, None, 0
    for a, b in itertools.permutations(reads, 2):
        length = overlap(a, b, min_length=k)
        if length > best_len:
            best_a, best_b, best_len = a, b, length
    return best_a, best_b, best_len

def greedy_scss(reads, k):
    while True:
        a, b, olen = pick_max_overlap(reads, k)
        if olen == 0:
            break
        reads.remove(a)
        reads.remove(b)
        reads.append(a + b[olen:])
    return ''.join(reads) # append all non-overlaps onto eachother and return the concatenated string
```

* **Question 3:** How many As are there in the full, assembled genome?

* **Question 4:** How many Ts are there in the full, assembled genome from the previous question?

Hint: the virus genome you are assembling is exactly 15,894 bases long

* Note: I chose the overlap length k=30 arbitrarily based on the range 30..50 mentioned during the lecture
for typical values used in the real-world.
```python
reads, _ = readFAST_Q('mystery.fq')
genome = greedy_scss(reads, 30)
print(len(genome))
print(genome.count('A'))
print(genome.count('T'))
```

### Supplemental
The greedy shortest common superstring algorithm is a
[Monte Carlo](https://en.wikipedia.org/wiki/Monte_Carlo_algorithm) algorithm:
```python
#
# optimal solution derived from the naive scss algorithm (which is very slow!)
#
res, _ = scss(['ABCD', 'CDBC', 'BCDA'])
print(res)
print(len(res))
```

```python
#
# non-optimal solution (ie. greedy solution is a monte-carlo algorithm [sometimes its wrong!])
#
res = greedy_scss(['ABCD', 'CDBC', 'BCDA'], 1)
print(res)
print(len(res))
```

### Solution 4
* [.ipynb](4_week/4_assignment.ipynb)
* [.py](4_week/4_assignment.py)
* [.html](https://nbviewer.jupyter.org/github/claytonjwong/Algorithms-DNA-Sequencing/blob/master/4_week/4_assignment.ipynb)
