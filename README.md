# Algorithms for DNA Sequencing
![](docs/double_helix.png)
* [Johns Hopkins University](https://www.jhu.edu/)
* [langmead-lab.org](http://www.langmead-lab.org/)

## 3 Laws of Assembly
1. If the suffix of read A is similar to the prefix of read B then A and B might overlap in the genome
    * See [Week 3 - Lecture 7](3_week/docs/overlaps_and_coverage.pdf) for details
2. More coverage leads to more and longer overlaps
    * See [Week 3 - Lecture 7](3_week/docs/overlaps_and_coverage.pdf) for details
3. Repeats make assembly difficult
    * See [Week 4 - Lecture 3](4_week/docs/third_law.pdf) for details 

## [Week 1: DNA Sequencing, Strings and Matching](1_week)
### Lectures
1. [Introduction](1_week/docs/intro.pdf)
2. [Why Study This?](1_week/docs/why_study_this.pdf)
3. [DNA sequencing past and present](1_week/docs/DNA_seq_past_present.pdf)
4. [Genomes as strings, reads as substrings](1_week/docs/genomes_strings.pdf)
5. [String definitions and Python examples](1_week/docs/python_str_def_ex.pdf)
6. [How DNA gets copied](1_week/docs/DNA_copying.pdf)
7. [How second-generation sequencers work](1_week/docs/second_gen_parallel.pdf)
8. [Sequencing errors and base qualities](1_week/docs/seq_errors_base_qualities.pdf)
9. [Working with sequencing reads (FASTQ format)](1_week/docs/fastq_format.pdf)
10. [Sequencers give pieces to genomic puzzles](1_week/docs/pieces_fragmentary.pdf)
11. [Read alignment and why it's hard](1_week/docs/read_alignment_hard.pdf)
12. [Naive exact matching](1_week/docs/naive_exact_matching.pdf)

### Assignment
* [.ipynb](1_week/1_assignment.ipynb)
* [.py](1_week/1_assignment.py)
* [.html](https://nbviewer.jupyter.org/github/claytonjwong/Algorithms-DNA-Sequencing/blob/master/1_week/1_assignment.ipynb)

### Resources
* [ERR037900_1.first1000.fastq](1_week/ERR037900_1.first1000.fastq)
* [lambda_virus.fa](1_week/lambda_virus.fa)

## [Week 2: Preprocessing, Indexing, and Approximate Matching](2_week)
### Lectures
1. [Boyer-Moore: Basics](2_week/docs/boyer-moore.pdf)
2. [Boyer-Moore: Putting It All Together](2_week/docs/boyer-moore-together.pdf)
3. [Diversion: Repetitive Elements](2_week/docs/repetitive_elements.pdf)
4. [Preprocessing](2_week/docs/preprocessing.pdf)
5. [Indexing and K-mers](2_week/docs/indexing_kmers.pdf)
6. [Data Structures for Indexing](2_week/docs/data_structures.pdf)
7. [Hash Tables for Indexing](2_week/docs/hash_tables.pdf)
8. [Variations on K-mer Indexes](2_week/docs/indexing_variations.pdf)
9. [Indexing by Suffix](2_week/docs/suffix.pdf)
10. [Approximate Matching, Hamming, and Edit Distance](2_week/docs/approximate.pdf)
11. [Pigeonhole Principle](2_week/docs/pigeonhole.pdf)

### Assignment
* [.ipynb](2_week/2_assignment.ipynb)
* [.py](2_week/2_assignment.py)
* [.html](https://nbviewer.jupyter.org/github/claytonjwong/Algorithms-DNA-Sequencing/blob/master/2_week/2_assignment.ipynb)

### Resources
* [bm_preproc.py](bm_preproc.py)
* [k-mer_index.py](k-mer_index.py)
* [chr1.GRCh38.excerpt.fasta](chr1.GRCh38.excerpt.fasta)

## [Week 3: Edit Distance, Assembly, and Overlaps](3_week)
### Lectures
1. [Edit Distance (part 1)](3_week/docs/edit_dist1.pdf)
2. [Edit Distance (part 2)](3_week/docs/edit_dist2.pdf)
3. [Edit Distance (part 3)](3_week/docs/edit_dist3.pdf)
4. [Edit Distance (part 4)](3_week/docs/edit_dist4.pdf)
5. [Global and Local Alignment](3_week/docs/global_and_local_alignment.pdf)
6. [De Novo Shotgun Assembly](3_week/docs/assembly_basics.pdf)
7. [1st and 2nd Laws of Assembly: Overlaps and Coverage](3_week/docs/overlaps_and_coverage.pdf)
8. [Overlap Graph](3_week/docs/overlap_graph.pdf)

### Assignment
* [.ipynb](3_week/3_assignment.ipynb)
* [.py](3_week/3_assignment.py)
* [.html](https://nbviewer.jupyter.org/github/claytonjwong/Algorithms-DNA-Sequencing/blob/master/3_week/3_assignment.ipynb)

### Resources
* [chr1.GRCh38.excerpt.fasta](chr1.GRCh38.excerpt.fasta)
* [ERR266411_1.for_asm.fastq](ERR266411_1.for_asm.fastq)

## [Week 4: Algorithms for Assembly](4_week)
### Lectures
1. [Shortest Common Substring](4_week/docs/scss.pdf)
2. [Greedy Shortest Common Substring](4_week/docs/greedy_scss.pdf)
3. [3rd Law of Assembly: Repeats are Bad](4_week/docs/third_law.pdf)
4. [De Bruijn Graphs and Eulerian Walks](4_week/docs/dbg1.pdf)
5. [When Eulerian Walks Go Wrong](4_week/docs/dbg2.pdf)
6. [Assemblers in Practice](4_week/docs/assemblers_in_practice.pdf)
7. [The Future is Long?](4_week/docs/longreads.pdf)
8. [Wrap Up](4_week/docs/wrap_up.pdf)

### Resources
* [mystery.fq](mystery.fq)

### Assignment
* [.ipynb](4_week/4_assignment.ipynb)
* [.py](4_week/4_assignment.py)
* [.html](https://nbviewer.jupyter.org/github/claytonjwong/Algorithms-DNA-Sequencing/blob/master/4_week/4_assignment.ipynb)

## Utility Functions
```python
def read_FAST_A(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
```

```python
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
```

## External Resources
* [Lectures](https://github.com/BenLangmead/ads1-slides)
* [Jupyter Notebooks](https://github.com/BenLangmead/ads1-notebooks)
* [Jupyter.org](https://jupyter.org/)
* [Fast Algorithms for Signal Processing](docs/Fast_Algorithms_for_Signal_Processing.pdf)

## Supplemental

* Jupyter Notebooks can be executed from the command line:

```
$ jupyter notebook 1_notebook.ipynb
[I 11:45:05.991 NotebookApp] The Jupyter Notebook is running at:
[I 11:45:05.991 NotebookApp] http://localhost:8889/?token=070644d6de70204df12235b2356476b577d0744b5df41422
[I 11:45:05.991 NotebookApp]  or http://127.0.0.1:8889/?token=070644d6de70204df12235b2356476b577d0744b5df41422
[I 11:45:05.991 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 11:45:05.998 NotebookApp] 
    
    To access the notebook, open this file in a browser:
        file:///Users/.../Library/Jupyter/runtime/nbserver-38748-open.html
    Or copy and paste one of these URLs:
        http://localhost:8889/?token=070644d6de70204df12235b2356476b577d0744b5df41422
     or http://127.0.0.1:8889/?token=070644d6de70204df12235b2356476b577d0744b5df41422
```

* The Python file for each Jupyter Notebook can be executed using ```ipython```.  If ```python``` is used
to execute then the following error will occur:
 
```
NameError: name 'get_ipython' is not defined
```

## Auxillary
K-mers are a fundamental concept for creating "words" from a DNA sequencing read.
These "words" are abstracted to computer science string algorithms (ie. simply finding pattern in text).

For example, a DNA substring consisting of two neucleotides is a 2-mer (regardless of Mr. Schwarzenegger's beliefs):

![](docs/its_not_a_2mer.png)
