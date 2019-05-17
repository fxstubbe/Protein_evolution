
### 1) load packages ###
# -------------------- #

library(data.table)
library(plyr)

# -------------------- #

### 2) path(s) and inputs ###
# ---------------------- #

f.root <- setwd("/Users/stubbf02/FX_Stubbe/ressources/proteins/staph")
f.gff <- paste( f.root , "/LukE/gff" , sep = "" )
f.save <- paste(f.root, "/LukE", sep = "")
setwd(f.gff)

gff.f <- "all_luke_aureus.gff3"

# -------------------- #

### 2) Fonctions ###
# ---------------------- #

gffRead <- function( gffFile , nrows = -1 ) {
  cat( "Reading " , gffFile , ": " , sep = "" )
  
  gff <- read.table( gffFile , sep = "\t" , as.is = TRUE , quote = "" ,
                     header = FALSE, comment.char = "#" , nrows = nrows ,
                     colClasses = c( "character" , "character" , "character" , "integer",  
                                     "integer", "character" , "character" , "character" , "character" ) )
  
  colnames( gff ) = c( "seqname" , "source" , "feature" , "start" , "end" ,
                       "score" , "strand" , "frame" , "attributes" )
  
  cat( "found" , nrow( gff ) , "rows with classes:" ,
       paste( sapply( gff , class ) , collapse = ", ") , "\n")
  stopifnot( !any( is.na( gff$start ) ) , !any( is.na( gff$end ) ) )
  return( gff )
}
getAttributeField <- function( x , field , attrsep = ";"){
  s = strsplit( x , split = attrsep , fixed = TRUE )
  sapply( s , function( atts ){
    a = strsplit( atts , split = "=", fixed = TRUE)
    m = match( field , sapply( a , "[" , 1 ))
    if ( !is.na( m ) ){
      rv = a[[ m ]][ 2 ]
    }
    else {
      rv = as.character( NA )
    }
    return( rv )
  } )
}


# -------------------- #

### 2) Cleaning dataframe ###
# ---------------------- #


gff <- gffRead( gffFile = gff.f )
gff <- gff[ -c( 3 , 4 , 5 , 6 , 7 , 8 ) ]


# remove lines of attributes that do not contain a strain ID
k <- 1
for( i in 1 : nrow(gff)) {
  
  if( grepl("strain", gff$attributes[k]) == FALSE) {
    gff <- gff[-k, , drop = FALSE]
  } else { k <- k + 1 }
}

# remove lines of source that is not genbank
k <- 1
for( i in 1 : nrow(gff)) {
  
  if( grepl("Genbank", gff$source[k]) == FALSE) {
    gff <- gff[-k, , drop = FALSE]
  } else { k <- k + 1 }
}

# extract attributes
gff$GenBank_locus_tag <- getAttributeField( gff$attributes , "strain")


# removing attributes
gff <- gff[ -c( 2 , 3 ) ]


# -------------------- #

### 3) writting output ###
# ---------------------- #

write.csv( gff, file = paste("gff_table_", gff.f,".csv", sep="") , row.names = FALSE , quote = FALSE) 



