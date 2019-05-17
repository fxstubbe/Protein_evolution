rm(list=ls())


library(stringr)

ReFF_name=c("USA300_FPR3757.fasta","USA500_2395.fasta","CFSAN007894.fasta","JH1.fasta")
ReFF_namelabels = sapply(ReFF_name, function(s) {
  s = gsub(".fasta", "", s)#replace ".fasta" with space 
})

genome.tested <-grep("300",ReFF_name)# find the one match "FPR"

#Path

f.root <- setwd('~/')
f.phylogeny <- paste(f.root, "/Users/stubbf02/Fx_Stubbe/ressources/proteins/staph/LukE", sep = "")
f.reference <- paste(f.root, "/FX_Stubbe/ressources/Prots/LukE/fasta", sep = "")
phylogeny_name = "SNP_matrix_3091_prots_95_cut_off.csv"
phyml_name = "SNP_matrix_.phy"
date = "03042019"

setwd( "/Users/stubbf02/Fx_Stubbe/ressources/proteins/staph/LukE" )
variants = read.csv(phylogeny_name, as.is = TRUE)

na.col <- which(as.numeric(sapply(variants, function(x) sum(is.na(x)))) >0)


# Assumes a fasta file representing a single genome
readFastaRef = function(refFile) {
  row = scan(refFile,what=character(0),sep="\n")
  chars = substr(row,1,1)
  base = chars!=">"
  seq = paste(row[base],collapse="")
  return(toupper(unlist(strsplit(seq,""))))
}
# Read the reference genome
setwd("/Users/stubbf02/Fx_Stubbe/ressources/proteins/staph/LukE/fasta/reference")
REF_Tested = readFastaRef("Usa300_FPR3757.fasta")

# Number of AA sites 
table(REF_Tested) 

# Identify the first column containing sequence data for the samples
#column 3 !!!

# Get the genome names
seqLabels = colnames(variants)[3:ncol(variants)]

# Count the number of prots
nSeq = length(seqLabels)

# Extract base calls & transpose so each row is a genome 
vSites = t(variants[,3:ncol(variants)])

# Ensure all characters are in upper case
vSites = toupper(vSites)

# Replace ambiguous calls
vSites[!(vSites=="A" | vSites=="R" | vSites=="D" | vSites=="N" | vSites=="C" | vSites=="E" | vSites=="Q" | vSites=="G" | vSites=="H" | vSites=="I" | vSites=="L" | vSites=="K" | vSites=="M" | vSites=="F" | vSites=="P" | vSites=="S" | vSites=="T" | vSites=="W" | vSites=="Y" | vSites=="V")] = "X"
#variants$Reference..base[!(variants$Reference..base=="A" | variants$Reference..base=="C" | variants$Reference..base=="G" | variants$Reference..base=="T")] = "N"


# Get the total length of the reference genome
lengthRef = length(REF_Tested)
lengthSites = length(vSites[1,])

# Isolate SNPs coordinates
vPos = sapply(variants[1], function(s) {
  as.numeric(unlist(s))
})
typeof(vPos)

# Sanity check
table(variants[,2],REF_Tested[vPos])


# Preparation of nucleotide sequence data into the correct format
# Specify PHYLIP file name
setwd("/Users/stubbf02/FX_Stubbe/ressources/proteins/staph/LukE")
phylipOutfile = paste(date,phyml_name, sep="")

# Output PHYLIP header: sequence number and sequence length
cat(paste(nSeq,lengthRef),sep="\n",file=phylipOutfile)

# For each genome, append to the PHYLIP file
#
for(i in 1:nSeq) {
  fullLengthRef = REF_Tested
  fullLengthRef[vPos] = vSites[i,]
  fullLengthRefCat = paste0(fullLengthRef,collapse="")
  cat(paste(seqLabels[i]," ",fullLengthRefCat,sep=""),sep="\n",file=phylipOutfile,append=TRUE)
  cat("Output",i,"of",nSeq,"\n")
}

paste("PhyML-3.1_macOS-MountainLion -i ",phylipOutfile," -b 0 -d aa -v 0 -c 1 -s BEST --no_memory_check", sep = "")
getwd()

##### On terminal : 
# PATH=$PATH:/Users/stubbf02/FX_Stubbe/ressources/tools/PhyML-3.1/
# PhyML-3.1_macOS-MountainLion -i 180125_annnotation_all_ST8.phylip -b 0 -v 0 -c 1 -s BEST --no_memory_check
# With bootstraps: PhyML-3.1_macOS-MountainLion -i 161102_phyML_ha_ac.phylip -b 1000 -v 0 -c 1 -s BEST --no_memory_check 