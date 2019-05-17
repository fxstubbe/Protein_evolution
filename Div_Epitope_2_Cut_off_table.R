##-------------------------------###
# Filter with cut-off
# Output table 
###-------------------------------###

### 0) Muscle align sequences ###
#In terminal
# -----------------------------------------------------------------------------------------------------------------------------#

muscle -in seqs.fa -out seqs_muscle.fa -maxiters 2


### 1) Load Packages ###
# -----------------------------------------------------------------------------------------------------------------------------#

library(stringr)
library(stringdist)
library(data.table)

### 2) Motif & Reference  ###
# -----------------------------------------------------------------------------------------------------------------------------#

ref.strain <- as.character(c("USA300_FPR3757"))
ref.motif <- as.character(c("YLPKNKIETTD"))

### 3) Path  ###
# -----------------------------------------------------------------------------------------------------------------------------#

f.root <- setwd( "/Users/stubbf02/FX_Stubbe/ressources/proteins/staph")
f.fasta <- paste( f.root , "/LukE/fasta" , sep = "" )
f.gff <- paste( f.root , "/LukE/gff" , sep = "" )
f.save <- paste( f.root , "/LukE" , sep = "" )

### 4) inputs  ###
# -----------------------------------------------------------------------------------------------------------------------------#

#Getting Gene tag & Strain (We removed strains that do not have the genbank attribute)
setwd( f.gff )
data.gff <- read.csv(file="gff_table_all_luke_aureus.gff3.csv", header=TRUE, sep=",")

#Importing the alignment
setwd( f.fasta )

row = scan("95_cutof_LukE_muscle.fa",what=character(0),sep="\n")
chars = substr(row,1,1)
seq_name <- row[grep(">",row)]

#Get the genban_ID, it's gonna be usefull to merge dataframe
Genbank_ID <- sapply( c( 1:length( seq_name ) ), function(k) {
  temp <- unlist( strsplit( seq_name[k], " " ) )[1]
  temp <- gsub( '>' , "", temp )
  return(temp)
  } 
)

#Get the sequences
sequences <- sapply( grep(">",row) , function(i) {paste(row[(i+1):(i+6)],collapse="")} )

#Storing Genbank_ID & sequewnces into a dataframe
data.fasta <-data.frame("seqname" = as.character(Genbank_ID),"sequences" = as.character(sequences))

#Merge data.gff & data.fasta by seqname
data.total <- merge(data.gff, data.fasta , by = "seqname")
names(data.total) <- c("Genbank_access", "strain", "sequences")


### 5) RMetrics identity between whole sequences ###
#Needs the reference prot sequence
# -----------------------------------------------------------------------------------------------------------------------------#

#Gets the row in dataframe where reference is
reference.protein <-grep(ref.strain ,data.total$strain)

#Add column containing levenshtein distance between prots
data.total$Prot_lv.distance <- sapply(data.total$sequences, function(k){
  stringdist(data.total$sequence[reference.protein] , k , method = "lv")
})

#Add column with %identity between prots
data.total$Prot_identity <- sapply(data.total$Prot_lv.distance, function(k){
  100 - ( ( k / nchar( as.character(data.total$sequence[reference.protein]))  * 100 ) )
})

### 6) Cut-Off ###
#Remove proteins that do not have at least 95% ( we only keep 50% of the strains like that !) of identity with the reference
#Why? Because there is a huge missanotation bias in the database
# -----------------------------------------------------------------------------------------------------------------------------#

cut_off = 95

data.total <- subset(data.total , Prot_identity >= cut_off, select=c(Genbank_access, strain, sequences, Prot_lv.distance, Prot_identity))

### 6) Writte_cut_off_table ###
# -----------------------------------------------------------------------------------------------------------------------------#

setwd( f.save )

write.csv( data.total, file = paste("LukE_95_cutof_realigned",".csv", sep="") , row.names = FALSE , quote = FALSE) 
