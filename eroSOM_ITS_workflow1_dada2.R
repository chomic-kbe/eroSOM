library("dada2")


# prepare data information
path <- "/mnt/data/chomic/seqme2021/eroSOM/its"
raw_reads <- sort(list.files("/mnt/data/chomic/seqme2021/eroSOM/its/raw", 
                             pattern="*.fq", full.names = TRUE))
sample_names <- sapply(strsplit(basename(raw_reads), "[.]"), `[`, 1)

# inspect quality profiles
plotQualityProfile(raw_reads[1:25])
plotQualityProfile(raw_reads[26:50])
plotQualityProfile(raw_reads[51:60])

# discard reads with unambiguous bases and >2 expected errors
filt_ee2 <- file.path(path, "filt_ee2", paste0(sample_names, "_filt_ee2.fq.gz"))
filt_ee2_out <- filterAndTrim(raw_reads, filt_ee2, maxN = 0, maxEE = 2, multithread = T)

# construct ASV table
err_ee2 <- learnErrors(filt_ee2, multithread = T)
dada_ee2 <- dada(filt_ee2, err=err_ee2, multithread=TRUE)
seqtab_ee2 <- makeSequenceTable(dada_ee2, orderBy = "abundance")

# remove chimeras 
seqtab_nochim_ee2 <- removeBimeraDenovo(seqtab_ee2, method="consensus", multithread=TRUE, 
                                        verbose=TRUE)
# write ASVs FASTA file for further taxonomic assignment
write.fasta.dada<-function(dada2, file){
  seqs<-dada2::getSequences(dada2)
  hash<-paste0(">",sapply(seqs, openssl::sha1, USE.NAMES = F))
  write(c(rbind(hash, seqs)),file)
}

write.fasta.dada(seqtab_nochim_ee2, paste0(path,"asvs.fa"))

# export ASV table
colnames(seqtab_nochim_ee2) <- sapply(colnames(seqtab_nochim_ee2), openssl::sha1, USE.NAMES = F)
rownames(seqtab_nochim_ee2) <- paste0("s", gsub("_filt_ee2.fq.gz", "", rownames(seqtab_nochim_ee2)))
head(seqtab_nochim_ee2)[1:6,1:6]

write.table(t(seqtab_nochim_ee2), paste0(path,"asv_table.txt"), quote = F)

