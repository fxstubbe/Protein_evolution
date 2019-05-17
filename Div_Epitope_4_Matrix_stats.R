
### 1) Load Packages ###
# ---------------------#

library( Biostrings )
library(stringr)
library(stringdist)
library(data.table)


### 2) Paths & Inputs  ###
# -----------------------#

f.root <- setwd( "/Users/stubbf02")
f.matrix <- paste( f.root , "/FX_Stubbe/ressources/Prots/LukE" , sep = "" )

matrix_name = "all_genbank_luke_aureus_matrix_no_filter.txt"
dataframe_name = "all_genbank_luke_aureus_dataframe.csv"

setwd( f.matrix )

variants = read.delim(matrix_name, as.is = TRUE)
data.gff <- read.csv(dataframe_name, header=TRUE, sep=",")

data.gff <- as.matrix(data.gff)
variants <- as.matrix(variants)


### 3)  Get the number of mutations compared to reference ###
# -----------------------#

table(data.gff[, 7])

### 3)  Remove lines with no mismatches###
# -----------------------#

k <- 1
for( i in 1 : nrow( variants ) ) {
  if( all(variants[k, 2] == variants [ k, 2 : 6209 ] ) == TRUE) {
    variants <- variants[-k, , drop = FALSE]
  } else { k <- k + 1 }
}

### 3)  Print table refering to mutations ###
# -----------------------#

for(i in 1:ncol(variants)){
  cat(paste("AA in reference : ", as.character(variants[i,2]), sep = ""))
  print(table(variants[i, ]))
}
  

write.table( variants, file = "all_genbank_luke_aureus_matrix.txt" , row.names = FALSE , quote = FALSE, sep = '\t') 

which(heatmap.protein$prot.occurence>4)
