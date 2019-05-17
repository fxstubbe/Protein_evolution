##-------------------------------###
# Makes a SNIP matrix from a MSA (muscle) file and a gff file containing name
###-------------------------------###

### 1) Muscle align sequences ###
#In terminal
# ---------------------#

muscle -in seqs.fa -out seqs_muscle.fa -maxiters 2


### 1) Load Packages ###
# ---------------------#

library( Biostrings )
library(stringr)
library(stringdist)
library(data.table)


### 2) Motif & Reference  ###
# -----------------------#

motif <- as.character(c("YLPKNKIETTD"))
ref <- as.character(c("USA300_FPR3757"))


### 3) Path  ###
# -----------------------#

f.root <- setwd( "/Users/stubbf02")
f.save <- paste( f.root , "/FX_Stubbe/ressources/Prots/LukE" , sep = "" )
f.fasta <- paste( f.root , "/FX_Stubbe/ressources/Prots/LukE/fasta" , sep = "" )
f.gff <- paste( f.root , "/FX_Stubbe/ressources/Prots/LukE/gff" , sep = "" )


### 4) inputs  ###
# -----------------------#

#Getting Gene tag & Strain 
setwd( f.gff )
data.gff <- read.csv(file="gff_table_all_luke_aureus.gff3.csv", header=TRUE, sep=",")

#Importing the alignment
setwd( f.fasta )

row = scan("allluke_aureus_muscle.fasta",what=character(0),sep="\n")
chars = substr(row,1,1)
seq_name <- row[grep(">",row)]

#Get the genban_ID, it's gonna be usefull to merge dataframe
Genbank_ID <- sapply(
  c(1:length(seq_name)), function(k){
    unlist(strsplit(seq_name[k], " "))[1]
  }
)

Genbank_ID <- sapply(
  c(1:length(seq_name)), function(k){
  str_remove(Genbank_ID[k], ">")})

#Get the sequences
sequences <- sapply( grep(">",row) , function(i) {paste(row[(i+1):(i+6)],collapse="")} )

#Storing Genbank_ID & sequewnces into a dataframe
data.fasta <-data.frame("seqname" = as.character(Genbank_ID),"sequences" = as.character(sequences))

#Merge data.gff & data.fasta by seqname
data.total <- merge(data.gff, data.fasta , by = "seqname")
names(data.total) <- c("Genbank_access", "strain", "sequences")


### 5) RMetrics identity between whole sequences ###
#Needs the reference prot sequence
# -----------------------#

#Gets the row in dataframe where reference is
reference.protein <-grep(ref ,data.total$strain)

#Add column containing levenshtein distance between prots
data.total$Prot_lv.distance <- sapply(data.total$sequences, function(k){
  stringdist(data.total$sequence[reference.protein] , k , method = "lv")
})

#Add column with %identity between prots
data.total$Prot_identity <- sapply(data.total$Prot_lv.distance, function(k){
  100 - ( ( k / nchar( as.character(data.total$sequence[reference.protein]))  * 100 ) )
})


### 6) Retrieving epitope motif ###
#SOME TWEAK HAS TO BE DONE HERE !!!!!

# -----------------------#
#Get vector with needed sequences

sequences <- unlist(as.character(data.total$sequence))

#Gets position of the motif

start <- str_locate(sequences[reference.protein], "LPKNKIETTD")[1]
start <- 148
end <- str_locate(sequences[reference.protein], "LPKNKIETTD")[2]

#Gets motifs
epitope <- sapply( c(1:nrow(data.total)) , function(k) { substr(unlist(as.character(data.total$sequences[k])), start, end) } )
data.total$epitope <- sapply( c( 1:length( epitope ) ), function( k ) { gsub('-', '', epitope[k] ) } ) 


### 7) ARMetrics identity between epitope motif  ###
# -----------------------#


#Add column containing levenshtein distance between motifs
data.total$motif_lv.distance <- sapply(data.total$epitope, function(k){
  stringdist(motif,k,method = "lv")
})

#Add column with %identity between motifs
data.total$motif_identity <- sapply(data.total$motif_lv.distance, function(k){
  100 - ( ( k / nchar( motif ) *100 ) )
})


### 8) SNP matrix ###
# -----------------------#

#Create a matrix where motifs are split (vewrtically)
#Add position in reference & reference amino_acid

SNP_matrix <- c()

SNP_matrix <- sapply(  1 : length( sequences )  , function(k) {
  unlist( strsplit( data.total$epitope[k] , "" ) )
})

SNP_matrix <- cbind( c(1 : 11) , unlist(strsplit(motif, "")) , SNP_matrix )


#Renaming columns with strain names
strain_names <- sapply(  1 : length( sequences )  , function(k) {
  as.character(data.total$strain[k]) 
})


colnames(SNP_matrix) <- c( "Motif_position" , "Ref_motif" , strain_names ) 


### 9) Outout filew ###
# -----------------------#

setwd( f.save )

write.table( SNP_matrix, file = "all_genbank_luke_aureus_matrix_no_filter.txt" , row.names = FALSE , quote = FALSE, sep = '\t') 
write.csv( data.total , file = "all_genbank_luke_aureus_dataframe.csv" , row.names = FALSE , quote = FALSE) 

