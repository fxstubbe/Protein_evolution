##-------------------------------###
# Get Identity Matrix
# Filter with cut-off
# Output table 
###-------------------------------###

### 1) Load Packages ###
# -----------------------------------------------------------------------------------------------------------------------------#

library(data.table)
library(tidyverse)

### 2) Motif & Reference  ###
# -----------------------------------------------------------------------------------------------------------------------------#

ref.strain <- as.character(c("USA300_FPR3757"))
ref.motif <- as.character(c("YLPKNKIETTD"))

### 3) Path ###
# -----------------------------------------------------------------------------------------------------------------------------#

f.root <- setwd( "/Users/stubbf02/FX_Stubbe/ressources/proteins/staph")
f.table <- paste( f.root , "/LukE/" , sep = "" )

### 4) Inputs ###
# -----------------------------------------------------------------------------------------------------------------------------#
setwd(f.table)

table <- fread( "LukE_95_cutof_realigned.csv" , sep = ",", header = TRUE , data.table = FALSE, fill = TRUE )
table$feature <- "LukE"


### 5) SNP matrix ###
# -----------------------#

#Gets the row in dataframe where reference is
reference.protein <-grep(ref.strain , table$strain)

#Create a matrix where motifs are split (vewrtically)
#Add position in reference & reference amino_acid

SNP_matrix <- c()

SNP_matrix <- sapply(  c(1 : nrow(table))   , function(k) {
  unlist( strsplit( table$sequences[k] , "" ) )
})

SNP_matrix <- cbind( c(1 : 311) ,  SNP_matrix )


#Renaming columns with strain names
strain_names <- sapply(  c(1 : nrow(table))  , function(k) {
  as.character(table$strain[k]) 
})
colnames(SNP_matrix) <- c( "Position", strain_names ) 

#Re ordering matrix to have reference in second column
SNP_matrix <- SNP_matrix[,c(1,3,2,4:ncol(SNP_matrix))]

#Remove lines with no mismatches (So you only keep lines where there is a SNP)
k <- 1
for( i in 1 : nrow( SNP_matrix ) ) {
  if( all(SNP_matrix[k, 2] == SNP_matrix [ k, 2 : ncol(SNP_matrix) ] ) == TRUE) {
    SNP_matrix <- SNP_matrix[-k, , drop = FALSE]
  } else { k <- k + 1 }
}

SNP_matrix <- data.frame(SNP_matrix)
dim(SNP_matrix)


### 5) Remove duplicate columns ###
# The dataset is big so making a phylogeny is gonna take forever
# Get unique sequences 
# -----------------------#

#Get column name
strain_names <- colnames(SNP_matrix)

#Collapse column into rows
LukE.sequences <- sapply( c( 2:3091 ), function( x ) { paste( SNP_matrix[,x], collapse = '' ) } )

#Get Plasmid positiom

Plasmid.Postition <- SNP_matrix$Position

#Make dataframe
SNP_matrix <- data.frame(cbind(strain_names[2:length(strain_names)], LukE.sequences))

#Get all_unique 
cleaned_SNP_matrix <- SNP_matrix %>% distinct(LukE.sequences, .keep_all = TRUE)

#Split row into columns (1 column by row, 1 row by character )
cleaned_sequences <- data.frame( sapply( c(1 : nrow( cleaned_SNP_matrix ) )   , function( k ) { unlist( strsplit( as.character( cleaned_SNP_matrix$LukE.sequences[k] ) , "" ) ) } ) )
colnames(cleaned_sequences) <- cleaned_SNP_matrix$V1
cleaned_SNP_matrix <- cbind( Plasmid.Postition , cleaned_sequences ) 


### 6) Output cleaned SNP_matrix ###
# -----------------------#

write.csv(cleaned_SNP_matrix, "SNP_matrix_3091_prots_95_cut_off_cleaned.csv")
