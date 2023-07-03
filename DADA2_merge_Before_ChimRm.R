#########################################################################################################################################

#title: DADA2 16S Processing Nudibranquios [Merge Oris Before Chimera Removal]
#author: Samuel Piquer-Esteban"
#date: 21/02/2022

#########################################################################################################################################

#0. Getting ready
print("#0. Getting ready")
print("Setting paths and parameters...")

# Set various files paths and changeable variables

##Set various paths

##A) Set paths to ori1 and ori2 results directories(input)
path_ori1="/storage/tbc/nudibranquios/DADA2/R_Objects_ori1"
path_ori2="/storage/tbc/nudibranquios/DADA2/R_Objects_ori2"

##B) Set path to Merge_Before_ChimRm (output)
path_out="/storage/tbc/nudibranquios/DADA2/R_Objects_merged/Merge_Before_ChimRm"

##Set changeable parameters

##A) Number of threadS
thrs=8

#Load Packages
print("Load packages ...")
library(dada2)
library(Biostrings)
#Done message
print("#0. Getting ready[Done]")

#########################################################################################################################################

#1.Merging Ori1 and Ori2 ASVs
print("#1. Merging Ori1 and Ori2 ASVs")

#Read tables
load(paste(path_ori1,"/seqtab_object.RData",sep=""))
seqtab_ori1=seqtab
load(paste(path_ori2,"/seqtab_object.RData",sep=""))
seqtab_ori2=seqtab

#Get ReverseComplement for Ori2
##RC Seqs function
rc_seqs<- function(i){
  require(Biostrings)
  rc_chr_string=as.character(reverseComplement(DNAString(i)))
  return(rc_chr_string)
}
##Get RC Sequences
seq_ori2_rc_seqs=sapply(colnames(seqtab_ori2), rc_seqs, USE.NAMES =F)
##Get copy of the matrix
seqtab_ori2_rc=seqtab_ori2
##Change seqs for RC Orientation
colnames(seqtab_ori2_rc)=seq_ori2_rc_seqs

#Merge tables 
seqtab_merged=mergeSequenceTables(table1 = seqtab_ori1, table2 = seqtab_ori2_rc, repeats = "sum")

#Print messages
print("seqtab_ori1 before merging:")
dim(seqtab_ori1)
print("seqtab_ori2 before merging:")
dim(seqtab_ori2)
print("seqtab_merged:")
dim(seqtab_merged)
# Save R objects
save(seqtab_merged,file=paste(path_out,"/seqtab_merged.RData",sep=""))

#Done message
print("Files created...")
print("#1. Merging Ori1 and Ori2 ASVs[Done]")

#########################################################################################################################################

#2. Identification and Removal of Chimeras
print("#2. Identification and Removal of Chimeras")

#Remove chimeras
seqtab_merged.nochim <- removeBimeraDenovo(seqtab_merged, method="consensus", multithread=thrs, verbose=TRUE)
# Save R objects
save(seqtab_merged.nochim,file=paste(path_out,"/seqtab_merged.nochim_object.RData",sep=""))

#Show info
print("Dimensions:")
dim(seqtab_merged.nochim)
print("Inspect distribution of sequence lengths:")
summary(nchar(getSequences(seqtab_merged.nochim)))
print("Percentage of Chimeric ASVs")
(1-(dim(seqtab_merged.nochim)[2]/dim(seqtab_merged)[2]))*100
print("Percentage of Chimeric count sequences")
(1-sum(seqtab_merged.nochim)/sum(seqtab_merged))*100
print("With respect to the previous step")

#Done message
print("Files created...")
print("#2. Identification and Removal of Chimeras[Done]")

#########################################################################################################################################

#3. Get Fasta and Count-table files
print("#3. Get Fasta and Count-table files")

# Giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab_merged.nochim)
asv_headers <- vector(dim(seqtab_merged.nochim)[2], mode="character")
for (i in 1:dim(seqtab_merged.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
# Making and writing out a fasta of our final ASV seqs
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste(path_out,"/ASVs_sequences.fa",sep=""))

# Writing Count table
asv_tab <- t(seqtab_merged.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, paste(path_out,"/ASVs_counts.tsv",sep=""), sep="\t", quote=F, col.names=NA)

#Done message
print("Files created...")
print("#3. Get Fasta and Count-table files[Done]")

#########################################################################################################################################

#4. SessionInfo
print("#4. SessionInfo")

sessionInfo()

print("#4. SessionInfo[Done]")
#Final message
print("Pipeline Done. Bye Dear...")

#########################################################################################################################################
