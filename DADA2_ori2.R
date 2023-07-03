#########################################################################################################################################

#title: DADA2 16S Processing ori2 Nudibranquios
#author: Samuel Piquer-Esteban"
#date: 17/02/2022

#########################################################################################################################################

#0. Getting ready
print("#0. Getting ready")
print("Setting paths and parameters...")

# Set various files paths and changeable variables

##Set various paths

##A) Set path to input FASTQ files (with-out PRIMERS and post filterbyname)
path1_R1 <- "/storage/tbc/nudibranquios/FilterbyName/R1_ori2"
path1_R2 <- "/storage/tbc/nudibranquios/FilterbyName/R2_ori2"

##Read in names for Forward and Reverse FASTQ filenames
fnFs <- sort(list.files(path1_R1, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path1_R2, pattern="_R2.fastq.gz", full.names = TRUE))
##Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_ori2_trimPrim_uniq"), `[`, 1)

##B) Set path to output DADA2 Clean FASTQ files
pathO1_R1="/storage/tbc/nudibranquios/clean_DADA2_fastqs/R1_ori2"
pathO1_R2="/storage/tbc/nudibranquios/clean_DADA2_fastqs/R2_ori2"
## Assign filenames for the filtered fastq.gz files
filtFs <- file.path(pathO1_R1, paste0(sample.names, "_ori2_trimPrim_uniq_DD2f_R1.fastq.gz"))
filtRs <- file.path(pathO1_R2, paste0(sample.names, "_ori2_trimPrim_uniq_DD2f_R2.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

##C) Set path to output R_Objects directory
pathO2="/storage/tbc/nudibranquios/DADA2/R_Objects_ori2"

##Set changeable parameters

##E) Number of threadS
thrs=8

##F) filterAndTrim() Changeable parameters
maxEE_R1=2
maxEE_R2=5
truncLen_R1=260
truncLen_R2=220
comp=TRUE

#Load Packages
print("Load packages ...")
library(dada2)
library(DECIPHER)
library(ShortRead)
library(Biostrings)

#Done message
print("#0. Getting ready[Done]")

#########################################################################################################################################

#1. Inspect read quality profiles
print("#1. Inspect read quality profiles")

# Self-made function
qc2pdf <- function(obj,out_name_path){
  #Get ranges
  ran=c(seq(1,length(obj),6),length(obj))
  #Open pdf
  pdf(out_name_path)
  #Loop get plots
  for(i in 1:(length(ran)-1)){
    if (i==(length(ran)-1)){
      a=ran[i]
      b=ran[i+1]
    } else {
      a=ran[i]
      b=ran[i+1]-1
    }
    #Save plot in pdf page
    print(plotQualityProfile(obj[a:b]))
  }
  #Save pdf file
  dev.off()
}
# Get quality profile plots Pre-QC with DD2 and save as PDF
qc2pdf(fnFs, paste(pathO2,"/ggplot_PreQC_R1.pdf",sep=""))
qc2pdf(fnRs, paste(pathO2,"/ggplot_PreQC_R2.pdf",sep=""))

#Done message
print("Files created...")
print("#1. Inspect read quality profiles[Done]")

#########################################################################################################################################

#2. Filter and quality trim
print("#2. Filter and quality trim")

# Process reads with DADA2
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(truncLen_R1,truncLen_R2),maxN=0, maxEE=c(maxEE_R1,maxEE_R2), truncQ=2, rm.phix=TRUE, compress=comp, multithread=thrs)
# Get quality profile plots Post-QC with DD2 ans save as PDF
qc2pdf(filtFs, paste(pathO2,"/ggplot_PostQC_R1.pdf",sep=""))
qc2pdf(filtRs, paste(pathO2,"/ggplot_PostQC_R2.pdf",sep=""))
#Save R object
save(out,file=paste(pathO2,"/filterAndTrim_out.RData",sep=""))

#Show filterAndTrim() Changeable parameters
#Show Dimensions
print("filterAndTrim() Changeable parameters:")
print(paste("maxEE_R1=",maxEE_R1))
print(paste("maxEE_R2=",maxEE_R2))
print(paste("truncLen_R1=",truncLen_R1))
print(paste("truncLen_R2=",truncLen_R2))

#Done message
print("Files created...")
print("#2. Filter and quality trim[Done]")

#########################################################################################################################################

#3. Generating error models
print("#3. Generating error models")

# Get error models for Forward and Reverse reads 
errF <- learnErrors(filtFs,multithread=thrs)
errR <- learnErrors(filtRs,multithread=thrs)
# Plot error models and save as PDF for Forward reads
pdf(paste(pathO2,"/ggplot_error_R1.pdf",sep=""))
print(plotErrors(errF, nominalQ=TRUE))
dev.off()
# Plot error models and save as PDF for Reads reads
pdf(paste(pathO2,"/ggplot_error_R2.pdf",sep=""))
print(plotErrors(errR, nominalQ=TRUE))
dev.off()
# Save R objects
save(errF,errR,file=paste(pathO2,"/error_models.RData",sep=""))

#Done message
print("Files created...")
print("#3. Generating error models[Done]")

#########################################################################################################################################

#4. Inferring ASVs
print("#4.Inferring ASVs")

# Inferring ASVs for Forward reads
dadaFs <- dada(filtFs, err=errF, pool=TRUE, multithread=thrs)
# Inferring ASVs for Reverse reads
dadaRs <- dada(filtRs, err=errR, pool=TRUE, multithread=thrs)
# Save R objects
save(dadaFs,dadaRs,file=paste(pathO2,"/dada_objects.RData",sep=""))

#Done message
print("Files created...")
print("#4.Inferring ASVs[Done]")

#########################################################################################################################################

#5. Merging forward and reverse reads
print("#5. Merging forward and reverse reads")

#Merge F and R sequences
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Save R objects
save(mergers,file=paste(pathO2,"/merge_object.RData",sep=""))

#Done message
print("Files created...")
print("#5. Merging forward and reverse reads[Done]")

#########################################################################################################################################

#6. Generating initial  Count-table
print("#6. Generating initial Count-table")

# Generate initial count table
seqtab <- makeSequenceTable(mergers)
# Save R objects
save(seqtab,file=paste(pathO2,"/seqtab_object.RData",sep=""))
#Show Dimensions
print("Dimensions:")
dim(seqtab)
# Inspect distribution of sequence lengths
print("Inspect distribution of sequence lengths:")
summary(nchar(getSequences(seqtab)))

#Done message
print("Files created...")
print("#6. Generating initial Count-table[Done]")

#########################################################################################################################################

#7. Identification and Removal of Chimeras
print("#7. Identification and Removal of Chimeras")

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=thrs, verbose=TRUE)
# Save R objects
save(seqtab.nochim,file=paste(pathO2,"/seqtabnochim_object.RData",sep=""))

#Show info
print("Dimensions:")
dim(seqtab.nochim)
print("Inspect distribution of sequence lengths:")
summary(nchar(getSequences(seqtab.nochim)))
print("Percentage of Chimeric ASVs")
(1-(dim(seqtab.nochim)[2]/dim(seqtab)[2]))*100
print("Percentage of Chimeric count sequences")
(1-sum(seqtab.nochim)/sum(seqtab))*100
print("With respect to the previous step")

#Done message
print("Files created...")
print("#7. Identification and Removal of Chimeras[Done]")

#########################################################################################################################################

#8. Get Fasta and Count-table files
print("#8. Get Fasta and Count-table Raw files")

# Giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
# Making and writing out a fasta of our final ASV seqs
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste(pathO2,"/ASVs_sequences_raw.fa",sep=""))

# Writing Count table
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, paste(pathO2,"/ASVs_counts_raw.tsv",sep=""), sep="\t", quote=F, col.names=NA)

#Done message
print("Files created...")
print("#8. Get Fasta and Count-table Raw files[Done]")

#########################################################################################################################################

#11. Get Track reads through the pipeline Table
print("#11. Track reads through the pipeline")

#Get pipeliline reads tracking as DADA2 Tutorial
## Note: If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input_Post_trimPrim_filterbyname", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

# Save R objects and Table
save(track,file=paste(pathO2,"/track_object.RData",sep=""))
write.table(track,paste(pathO2,"/track_table.tsv",sep=""),sep="\t",quote=F,col.names=NA)

#Done message
print("Summary files created...")
print("#11. Track reads through the pipeline[Done]")

#########################################################################################################################################

#12. SessionInfo
print("#12. SessionInfo")

sessionInfo()

print("#12. SessionInfo[Done]")
#Final message
print("Pipeline Done. Bye Dear...")

#########################################################################################################################################
