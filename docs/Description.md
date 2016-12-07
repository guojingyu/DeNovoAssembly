##Description
This project is a simple implementation for solving DNA sequence de novo assembly. It is implemented in python 2. 

###Implementation Assumptions



###Description of the project
This project is implemented with refering to kmer built De Bruijn Graph (DBG) 
and Euler Walk method, which is to find a path/circuit to traverse all edges in 
the graph.  

The DBG method avoided (to some extent) overcollapsed repeats problem that might be introduced by the overlapping based methods (overlap-layout-consensus). As well it is used more widely in NGS based assembly applications, thus it was chosen for this project.

####Kmer, De Bruijin Graph and Eulerian

Each kmer, as a substring of a DNA sequence of fixed length of k 
(positive integer), represents an edge in the DBG. The prefix and suffix 
(k-1)mer of each kmer form two nodes for the edge. 

DBG is formed as a multi directed graph, which allows more than one edges 
between a pair of nodes or more than one loop on a single nodes. Although 
permited multiple edges, in DBG, the nodes ((k-1)mers) are unique (in terms 
of DNA sequence/string it represents).

Eulerian graph is true if and only if the in degree and out degree of each 
nodes in the graph is equal. When such criteria is satisified, there exists 
(at least) one circuit in the graph that pass through each and every edge only 
once. This provides us a manner to build a circular sequence from 
the edges (or the nodes) in the order of traverse. However, it is possible 
that in a (strongly) connected directed graph, one node could have one extra out degree and 
another could have one extra in degree, while the rest nodes remains equal of
 in and out degree. These two unbalanced nodes may represents a start and an
  end for a Eulerian path, which would be, in this exercise, a linear 
  DNA sequence.

###Results
For the dummy_data.fasta file, a 19 bp long DNA assembly can be found with a k > 5, as **ATTAGACCTGCCGGAATAC**, which is the same as reference answer.
 
For the provided 50 fasta sample file, when k is set to 12-15, the output assembled DNA sequence has closely ranged around 20000 bp long but not unique. When the k is set to larger than 20, it would be more stablized to a sequence of **19914** bp 
long. 

It is also noticeable, for the provided 50 fasta sample file, that the 
running time increased when k is getting larger. This holds true at least 
when k is still relatively small (20-30) than the ~1000 length of the fasta 
records. It may take less than 1 minute to finish if k = 20, but 
more than 1 minute to finish at 30 (single threaded on a typical laptop with
 i7 processor). 

###Resources
1. Phillip E C Compeau, Pavel A Pevzner, Glenn Tesler. How to apply de Bruijn graphs to genome assembly. Nature Biotechnology 29, 987–991 (2011) doi:10.1038/nbt.2023
2. Zhenyu Li, Yanxiang Chen, Desheng Mu, Jianying Yuan, Yujian Shi, Hao Zhang, Jun Gan, Nan Li, Xuesong Hu, Binghang Liu, Bicheng Yang, and Wei Fan. Comparison of the two major classes of assembly algorithms: overlap–layout–consensus and de-bruijn-graph. Briefings in Functional Genomics. 2011 doi:10.1093/bfgp/elr035
3. https://www.coursera.org/learn/genome-sequencing/home/  (UCSD coursera 
Bioinformatics course)