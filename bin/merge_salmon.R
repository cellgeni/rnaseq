#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

tpms <- NULL
counts <- NULL
for(f in head(args, -1)) {
    d <- read.table(f, header = T)
    tpms <- cbind(tpms, d[,4])
    counts <- cbind(counts, d[,5])
}

colnames(tpms) <- head(args, -1)
colnames(counts) <- head(args, -1)

rownames(tpms) <- d[,1]
rownames(counts) <- d[,1]

write.csv(
    tpms, 
    file = paste0("merged_tpms_", tail(args, 1), ".csv"), 
    quote = F,
    row.names = F
)

write.csv(
    counts, 
    file = paste0("merged_counts_", tail(args, 1), ".csv"), 
    quote = F,
    row.names = F
)
