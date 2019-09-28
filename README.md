# Algorithms for DNA Sequencing

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
    * Unrelated humans have genomes of ~3 billion base pairs which are 99.8 - 99.9% similar
11. [Read alignment and why it's hard](1_week/docs/read_alignment_hard.pdf)
12. [Naive exact matching](1_week/docs/naive_exact_matching.pdf)

### Notes and Assignments
* [1_notebook.ipynb](1_week/1_notebook.ipynb)
    * [.md](1_week/1_notebook.md)
    * [.html](1_week/1_notebook.html)
    
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

### Notes and Assignments
* [2_notebook.ipynb](2_week/2_notebook.ipynb) 
    * [.md](2_week/2_notebook.md)
    * [.html](2_week/2_notebook.html)
    
### Resources
* [bm_preproc.py](bm_preproc.py)
* [k-mer_index.py](k-mer_index.py)
* [chr1.GRCh38.excerpt.fasta](chr1.GRCh38.excerpt.fasta)

## External Resources
* [Lectures](https://github.com/BenLangmead/ads1-slides)
* [Jupyter Notebooks](https://github.com/BenLangmead/ads1-notebooks)
* [Jupyter.org](https://jupyter.org/)

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