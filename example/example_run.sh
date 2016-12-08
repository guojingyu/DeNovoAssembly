#!/bin/sh

# assuming in the project folder DeNovoAssembly

# To run the dummy_data.fasta:
python ./de_novo_assembly/run.py -i ./data/dummy_data.fasta -k 6 -o ./output/dummy_data_assembly_output.txt --print_to_console

# To include a graph, simply add the option as '--graph':
python ./de_novo_assembly/run.py -i ./data/dummy_data.fasta -k 6 -o ./output/dummy_data_assembly_output.txt --print_to_console --graph

# To run the 50 fasta record set:
python ./de_novo_assembly/run.py -i ./data/coding_challenge_data_set.fasta -k 20 -o ./output/test_data_assembly_output.txt