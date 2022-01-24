library(seqinr)
library(LncFinder)

args <- commandArgs(TRUE)
list_dir <- list.dirs(path = args[1], recursive = FALSE)
for (w in list_dir) {

lnc <- read.fasta(args[2])
mrna <- read.fasta(paste(w, '/train_mrna.fasta', sep=''))
strc <- args[3]
print(strc)

if (strc == "SS") {
print('secondary structure on')
RNAfold.path <- '/home/pronozinau/miniconda3/bin/RNAfold'

mrna1 <- run_RNAfold(mrna, RNAfold.path = RNAfold.path, parallel.cores = 40)
lnc1 <- run_RNAfold(lnc, RNAfold.path = RNAfold.path, parallel.cores = 40)

feach <- make_frequencies(cds.seq = mrna1,  mRNA.seq = mrna1, lncRNA.seq = lnc1,
SS.features = TRUE, cds.format = "SS",
lnc.format = "SS", check.cds = TRUE,
ignore.illegal = FALSE)
print('done')
saveRDS(feach, file=paste(w, '/frequencies.rds', sep=""))

my_model <- build_model(lncRNA.seq = lnc1, mRNA.seq = mrna1,
frequencies.file = feach, SS.features = TRUE,
lncRNA.format = "SS", mRNA.format = "SS",
parallel.cores = 2, folds.num = 2, seed = 1)
saveRDS(my_model, file=paste(w,'/model.rds',sep=""))

print(paste(w, '/test.fasta', sep=''))
seq <- read.fasta(paste(w, '/test.fasta', sep=''))
test <- run_RNAfold(seq, RNAfold.path = RNAfold.path, parallel.cores = 40)
result_2 <- lnc_finder(test, SS.features = TRUE, format = "SS",
frequencies.file = feach, svm.model = my_model,
parallel.cores = 2)
write.table(result_2, file=paste(w, '/results.csv', sep=""), col.names=F, sep = ',')
} else {
print('secondary structure off')

feach <- make_frequencies(cds.seq = mrna,  mRNA.seq = mrna, lncRNA.seq = lnc,
SS.features = FALSE, cds.format = "DNA",
lnc.format = "DNA", check.cds = TRUE,
ignore.illegal = FALSE)
print('done')
saveRDS(feach, file=paste(w, '/frequencies.rds', sep=""))

my_model <- build_model(lncRNA.seq = lnc, mRNA.seq = mrna,
frequencies.file = feach, SS.features = FALSE,
lncRNA.format = "DNA", mRNA.format = "DNA",
parallel.cores = 2, folds.num = 2, seed = 1)
saveRDS(my_model, file=paste(w,'/model.rds',sep=""))

print(paste(w, '/test.fasta', sep=''))
seq <- read.fasta(paste(w, '/test.fasta', sep=''))
result_2 <- lnc_finder(seq, SS.features = FALSE, format = "DNA",
frequencies.file = feach, svm.model = my_model,
parallel.cores = 2)
write.table(result_2, file=paste(w, '/results.csv', sep=""), col.names=F, sep = ',')
}
}
print(list_dir)