#!/bin/bash

FIN=$1 # input phylip format distance matrix
FOUT=$2 # output newick format nj tree

## fastme is installed using anaconda

fastme -i $FIN -o $FOUT -n -s
