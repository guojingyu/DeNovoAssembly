## DNA Sequence De Novo Assembly

This repo contains an implementation of a solution to DNA sequence assembly problem. The is based on Euler Walk on De Bruijn Graph constructed using kmer representation. 

(For further details, please read docs/Description.md)

### Prerequisites
The following only represents one possible set of versions of required libraries, which is used in development.
```
python == 2.7.11
nose == 1.3.7
numpy == 1.11.1
biopython == 1.68
networkx == 1.11
matplotlib == 1.5.3
```

### How to run the code
To run the dummy_data.fasta:  
```
python ./de_novo_assembly/run.py -i ./data/dummy_data.fasta -k 5 -o ./output/dummy_data_assembly_output.txt --print_to_console
```

To include a graph, simply add the option as '--graph':  
```
python ./de_novo_assembly/run.py -i ./data/dummy_data.fasta -k 5 -o ./output/dummy_data_assembly_output.txt --print_to_console --graph
```

To run the 50 fasta record set:  
```
python ./de_novo_assembly/run.py -i ./data/coding_challenge_data_set.fasta -k 21 -o ./output/test_data_assembly_output.txt
```

The k is set to 5 for the dummy dataset (works when 4<k<=9), and 20 for the 50 fasta record dataset for stable result. The graph plotting feature is very primitive and may not work well for large De Bruijn Graph visualization. A sample plot for the dummy data of 4 fasta records can be found here for a 15 node 5mer De Bruijn Graph for a total 19 bp assembly.

![Image of dummy dataset DBG](https://github.com/guojingyu/DeNovoAssembly/blob/master/dummy_data_de_bruijn.png)

For all running options:  
```
DeNovoAssembly $ python ./de_novo_assembly/run.py -h
usage: run.py [-h] [-i INPUTFASTAFILE] [--reverse_complement] [-k KMERLENGTH]
              [--graph] [-o OUTPUTFILE] [--output_longest_assembly]
              [--print_to_console]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFASTAFILE, --inputfastafile INPUTFASTAFILE
                        fasta formatted text file containing the sequences for
                        assembly
  --reverse_complement  include or not. to include to the reverse complement
                        sequences of the fasta file sequence
  -k KMERLENGTH, --kmerlength KMERLENGTH
                        parameter k to control the length of kmer as the edge
                        of De Bruijn Graph
  --graph               include or not. plot the De Bruijn Graph constructed
                        from kmer of the input sequence
  -o OUTPUTFILE, --outputfile OUTPUTFILE
                        output assembled DNA sequences into file with input
                        name and time stamp appended
  --output_longest_assembly
                        include or not. output the longest assembled DNA
                        sequence. If more than one same longest DNA sequence
                        found, return the first one in rank
  --print_to_console    include or not. to print assembled DNA sequence to
                        console

```



