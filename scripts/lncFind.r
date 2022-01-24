library(seqinr)
library(LncFinder)

args <- commandArgs(TRUE)
path <- readLines(args[1])
strc <- args[4]
print(strc)
if (strc == "SS") {

frequencies <- readRDS(paste(path, 'frequencies.rds', sep=""))
model <- readRDS(paste(path, 'model.rds', sep=""))

if (file_test("-d", args[2]) == TRUE) {
file <- list.files(path = args[2], pattern = ".fa", full.names=TRUE)
for (w in file) {
print(w)
seq <- read.fasta(w)
test <- run_RNAfold(seq, RNAfold.path = RNAfold.path, parallel.cores = 40)
result_2 <- lnc_finder(test, SS.features = TRUE, format = "SS",
frequencies.file = frequencies, svm.model = model,
parallel.cores = 2)
write.table(result_2, file=paste(args[3],'/lncFinder_train.csv',sep=""), append=T, col.names=F, sep = ",")
print("loop done")
}
} else {
print(args[2])
seq <- read.fasta(args[2])
test <- run_RNAfold(seq, RNAfold.path = RNAfold.path, parallel.cores = 40)
result_2 <- lnc_finder(test, SS.features = TRUE, format = "SS",
frequencies.file = frequencies, svm.model = model,
parallel.cores = 2)
write.table(result_2, file=paste(args[3],'/lncFinder_train.csv',sep=""), col.names=F, sep = ',')
print("no loop")
}
} else {
frequencies <- readRDS(paste(path, 'frequencies.rds', sep=""))
model <- readRDS(paste(path, 'model.rds', sep=""))

if (file_test("-d", args[2]) == TRUE) {
file <- list.files(path = args[2], pattern = ".fa", full.names=TRUE)
for (w in file) {
print(w)
seq <- read.fasta(w)
result_2 <- lnc_finder(seq, SS.features = FALSE, format = "DNA",
frequencies.file = frequencies, svm.model = model,
parallel.cores = 2)
write.table(result_2, file=paste(args[3],'/lncFinder_train.csv',sep=""), append=T, col.names=F, sep = ",")
print("loop done")
}
} else {
print(args[2])
seq <- read.fasta(args[2])
result_2 <- lnc_finder(seq, SS.features = FALSE, format = "DNA",
frequencies.file = frequencies, svm.model = model,
parallel.cores = 2)
write.table(result_2, file=paste(args[3],'/lncFinder_train.csv',sep=""), col.names=F, sep = ',')
print("no loop")
}
}