#!/usr/bin/bash

random_gene=$(shuf -i 1-4100 -n 1)

# Run R implementation
echo "Run R implementation on gene ${random_gene}"
Rscript test/test.R ${random_gene}

# Run python implementation
# Since python starts counting from 0, we subtract 1 from the randomly
# selected gene
echo "Run Python implementation on gene ${random_gene}"
PYTHONPATH=./python python2 test/test.py $((${random_gene}-1))

# Check resulting computations
Rscript test/verify.R
