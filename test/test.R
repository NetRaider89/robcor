source('R/cormad.R')

argv = commandArgs(trailingOnly=TRUE)
selected_gene = as.integer(argv[1])

data = as.matrix(read.table('data/RNASeq.txt'))
n_samples = dim(data)[1]

correlations = list()
for(i in seq(n_samples)){
    correlations = c(correlations, cormad(data[selected_gene, ], data[i, ]))
}

write.table(correlations, 'test/R_correlations.txt', quote=FALSE, sep="\n", 
            row.names=FALSE, col.names=FALSE)
