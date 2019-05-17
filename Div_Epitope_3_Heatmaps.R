##-------------------------------###
# Get Identity Matrix
# Filter with cut-off
# Output table 
###-------------------------------###

### 1) Load Packages ###
# -----------------------------------------------------------------------------------------------------------------------------#

library(stringr)
library(stringdist)
library(data.table)
library(ggplot2)
library(RColorBrewer)

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

### 5) Make heatmap Matrix###
# -----------------------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------Get the reference strain

#Gets the row in dataframe where reference is
ref.protein <-table[grep(ref.strain ,table$strain),]

#Remove reference row in dataframe (because we compare to this particular protein)
table <- table[-grep(ref.strain ,table$strain),]

#----------------------------------------------------Occurence calculation for the complete protein

#Data Wrangling, make a table that will be used for the dataframe
ref.position <- c( 1:nchar( ref.protein$sequences ) )
prot.occurence <- sapply( c( 1:nchar( ref.protein$sequences ) ), function(x) {
  
  ref.position <- substr(ref.protein$sequences, x, x)
  occurence <- 0
  
  for ( i in 1:nrow( table ) ) {
    tested.position <- substr(table$sequences[i], x, x)
    if (ref.position  != tested.position){
      occurence <- occurence + 1
    }
  } 
  return ( occurence ) 
  } )


heatmap.protein <- data.frame( cbind( ref.position, prot.occurence  ) )
heatmap.protein$feature <- "LukE"
heatmap.protein <- heatmap.protein[,c(3,1:2)]

#----------------------------------------------------Occurence calculation for the epitope

#Get the epitope
start <- str_locate(ref.protein$sequence, "YLPKNKIETTD")[1]
end <- str_locate(ref.protein$sequence, "YLPKNKIETTD")[2]
ref.epitope <- substr(ref.protein$sequences, start, end)


#Get the epitope 
table$epitope_seq <- sapply( c( 1:nrow( table ) ), function(x){
  substr(table$sequences[x], start, end)
} ) 

#occurence calculation
ref.epitope.position <- c( 1:nchar( ref.epitope ) )
epitope.occurence <- sapply( c( 1:nchar( ref.epitope ) ), function(x) {
  
  ref.position <- substr(ref.epitope, x, x)
  occurence <- 0
  
  for ( i in 1:nrow( table ) ) {
    tested.position <- substr(table$epitope_seq[i], x, x)
    if (ref.position  != tested.position){
      occurence <- occurence + 1
    }
  } 
  return ( occurence ) 
} )

heatmap.epitope <- data.frame( cbind( ref.epitope.position, epitope.occurence  ) )
heatmap.epitope$feature <- "LukE"
heatmap.epitope <- heatmap.epitope[,c(3,1:2)]


### 5) Plot the heatmap  for the whole protein ###
# -----------------------------------------------------------------------------------------------------------------------------#

display.brewer.all()

for(i in 148:160){
  heatmap.protein[i,3] <- NA
}

#create a new variable from incidence
heatmap.protein$OccurenceFactor <- cut(heatmap.protein$prot.occurence,
                          breaks = c(-1,0,1,10,100),labels=c("0","0-1","1-10","10-100"))

#change level order
heatmap.protein$OccurenceFactor <- factor(as.character(heatmap.protein$OccurenceFactor),
                             levels=rev(levels(heatmap.protein$OccurenceFactor)))



#define a colour for fonts
textcol <- "grey40"

p <- ggplot(heatmap.protein,aes(x= ref.position,y=feature,fill=OccurenceFactor))+ 
  geom_tile(colour="white",size=0.01)+
  #remove axis labels, add title
  labs(x="",y="",title="") +
  #remove extra space
  scale_y_discrete(expand=c(0,0))+
  #custom breaks on x-axis
  #scale_x_discrete(expand=c(0,0),
                   #breaks=c(0,50,100,150,200,250,300,350))+
  #custom colours for cut levels and na values
  scale_fill_manual(expand=c(0,0),values=c("#abdda4","#fdae61","#fee08b","#d53e4f"),na.value="#3F7FBF")+
  theme_minimal(base_size = 10)

p

### 5) Plot the heatmap  for the epitope ###
# -----------------------------------------------------------------------------------------------------------------------------#

heatmap.epitope$OccurenceFactor <- cut(heatmap.epitope$epitope.occurence,
                                       breaks = c(-1,0,1,10,100),labels=c("0","0-1","1-10","10-100"))

#change level order
heatmap.epitope$OccurenceFactor <- factor(as.character(heatmap.epitope$OccurenceFactor),
                                        levels=rev(levels(heatmap.epitope$OccurenceFactor)))



#define a colour for fonts
textcol <- "grey40"

g <- ggplot(heatmap.epitope,aes(x= ref.epitope.position,y=feature,fill=OccurenceFactor))+ 
  geom_tile(colour="white",size=0.01)+
  #remove axis labels, add title
  labs(x="",y="",title="") +
  #remove extra space
  scale_y_discrete(expand=c(0,0))+
  #custom breaks on x-axis
  #scale_x_discrete(expand=c(0,0),
  #breaks=c(0,50,100,150,200,250,300,350))+
  #custom colours for cut levels and na values
  scale_fill_manual(expand=c(0,0),values=c("#abdda4","#abdda4","#abdda4","#abdda4"),na.value="grey90")+
  #scale_fill_manual(values=rev(brewer.pal(7,"OrRd")),na.value="grey90")+
  theme_minimal(base_size = 10)

g

which(heatmap.protein$prot.occurence>4)
