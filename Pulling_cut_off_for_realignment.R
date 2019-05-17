##-------------------------------###
# Extract sequnces to re-muscle-align
###-------------------------------###

### 1) Load Packages ###
# -----------------------------------------------------------------------------------------------------------------------------#

library(stringr)
library(stringdist)
library(data.table)

### 2) Path  ###
# -----------------------------------------------------------------------------------------------------------------------------#

f.root <- setwd( "/Users/stubbf02/FX_Stubbe/ressources/proteins/staph/LukE")
f.fasta <- paste( f.root , "/fasta" , sep = "" )
f.gff <- paste( f.root , "/gff" , sep = "" )


### 3) inputs  ###
# -----------------------------------------------------------------------------------------------------------------------------#

#Importing the FASTA to clean
setwd( f.fasta )
row = scan("all_luke_aureus.fasta",what=character(0),sep="\n")
chars = substr(row,1,1)
seq_name <- row[grep(">",row)]
sequences <- sapply( grep(">",row) , function(i) {paste(row[(i+1):(i+5)],collapse="")} )
table <- data.frame(cbind(seq_name, sequences))

#Importing the cut_off table
setwd( f.root )
data.cutof <- read.csv(file="LukE_95_cutof.csv", header=TRUE, sep=",")

### 4) Output  ###
# -----------------------------------------------------------------------------------------------------------------------------#

setwd( f.fasta )

output.file = paste( "95_cutof_LukE.fasta" , sep = '' )

for ( i in 1:nrow( data.cutof ) ){
  cat( paste( as.character(table[ grep(as.character(data.cutof$Genbank_access[i]), table$seq_name ) , 1 ] ), sep = "" ), sep = "\n", file = output.file, append = TRUE)
  cat( paste( as.character(table[ grep(as.character(data.cutof$Genbank_access[i]), table$seq_name ) , 2 ] ), sep = "" ), sep = "\n", file = output.file, append = TRUE)
}

### 0) Muscle align sequences ###
#In terminal
# -----------------------------------------------------------------------------------------------------------------------------#

cd /Users/stubbf02/FX_Stubbe/ressources/proteins/staph/LukE/fasta
muscle -in 95_cutof_LukE.fasta -out 95_cutof_LukE_muscle.fa 



nchar(as.character(table[ grep("FPR", table$seq_name ), 2]))
       