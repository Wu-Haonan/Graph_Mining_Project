# Graph_Mining_Project

Course Project of Graph Mining

This repo used for course project.

## Problem statement

**V(D)J recombination** is a fundamental genetic mechanism that enables the adaptive immune system to generate an extraordinarily diverse repertoire of antigen receptors on B and T lymphocytes. By randomly assembling variable (V), diversity (D), and joining (J) gene segments and undergoing somatic hypermutation, these processes produc a vast array of immunoglobulin (Ig) and T cell receptor (TCR) sequences, each with unique antigen-binding specificity. This molecular diversity underlies the immune system’s ability to recognize and respond to an enormous range of pathogens. Accurately reconstructing the immune repertoire—i.e., the full set of antigen receptor sequences present in an individual—is essential for understanding immune responses, tracking clonal expansions in infection or cancer, and developing targeted immunotherapies. High-throughput sequencing and computational reconstruction of V(D)J repertoires have therefore become indispensable tools in immunology and translational medicine.  



Previous methods depend on the known information of V, D, and J sequences; they first map sequences in the repertoire to the reference and then determine the class of V and J genes. Then, combining more biological knowledge to infer the clonal lineage of the repertoire. In this course project, I try to build a de novo pipeline, i.e., employing the sequence distance and calling a graph-clustering algorithm (we choose [Leidan](https://leidenalg.readthedocs.io/en/stable/intro.html) here). The details can be found in the following section. 



## Ground Truth

1. The original VDJ sequences are stored in `./data/PRJNA324093_Dnr4_10k.fasta`

2. ImmunoTools (https://immunotools.github.io/immunotools/) is used to align and analyze VDJ sequences. The results are stored in `./Alignment_results`. The output we focused on is `compressed_cdr3s.fasta`, which shows all the CDR3 sequences, counts, and corresponding V, J genes. 

3. Run `./clonal_lineage` to generate the ground truth of classification in the CSV file `ground_truth.csv`. Based on the following rules.

   a. they have the same CDR3 length,
   b. the Hamming distance between their nucleotide sequences is below 10%
   c. They represent VDJ sequences with identical V and J genes (ignoring gene alleles).

## Pipeline

### All-vs-all alignment

We call [EdLib](https://martinsos.github.io/edlib/) to implement all-vs-all global alignment. If two strings have edit distance $\leq 40$, the two nodes will have an edge and assign the distance as weight for the edge. 

### Leiden cluster

We call Leiden algorithms to implement clustering and output the results in `Leiden_cluster.csv`.

### Evaluate

We use ARI to evaluate the clustering results. The results are shown below

```tex
Ground truth cluster count: 403
Leiden cluster count: 151
Adjusted Rand Index (ARI): -0.0004
```



## Refection

The results are not good and close to random. I think it can be interpreted by the following reason

1. I set the threshold $40$ to connect an edge, which lead to a at least $90\%$​​ identity for most sequences. This threshold may too high, so we cannot capture complete information to cluster. 
1. The Leiden algorithm may not be suitable for this problem. 



For further improvements, we can extend it in the following directions.

1. Try to determine a better threshold to connect edges. 
2. Try to explore a faster way to speed up the all-vs-all distance calculations. 
3. Try to develop novel graph-clustering algorithms to adapt to the clonal lineage classification problem. 

