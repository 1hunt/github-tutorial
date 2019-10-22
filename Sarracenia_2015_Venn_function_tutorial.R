###2015 Sarracenia Data Venn Diagram Comparisons###

##Load relevant packages
if (!require("eulerr")) {install.packages("eulerr"); require("eulerr")}
if (!require("vegan")) {install.packages("vegan"); require("vegan")}

##Load Data
setwd("/Users/Jacob/Desktop/Bittleston/Sarracenia_Data/Vennpdfs/")
#Raw Data
sar <- read.table("asv-table-dada2-18s-2015-sarracenia.txt", header = T, row.names = 1, check.names = F, sep = "\t") 
sar <- sar[order(row.names(sar)),]
#Metadata
mdat <- read.csv("Metadata_2015_Sarracenia.csv", head=T, row.names=1) 
row.names(mdat) == row.names(sar)
mdat <- mdat[-38,] #sar(18sdata) does not include the row found here
sar <- subset(sar, mdat$Project %in% ("Comparison"))
mdat <- subset(mdat, row.names(mdat) %in% row.names(sar))
mdat <- mdat[order(row.names(mdat)),]
row.names(mdat) == row.names(sar)

###Function: sar_venn()
##Syntax: sar_venn(species,type,remove.n)
##Arguments: 
#species = any combination of "alata", "flava", "leucophylla", "purpurea", "rosea", "rubra"
#type = character, "euler" or "venn", default="euler"
#remove.n = numeric, default=0, removes x number of OTUs off the bottom

sar_venn <- function(species,type="euler",remove.n=0){
  
#Raw data
sar_sbst <- subset(sar, mdat$Host_Spp %in% species)
#Metadata
mdat_sbst <- subset(mdat, row.names(mdat) %in% row.names(sar_sbst))
#Drop old factors from metadata
mdat_sbst$Host_Spp <- mdat_sbst$Host_Spp[drop=T]

#Rarify
sar_sbst.r <- rrarefy(sar_sbst, min(rowSums(sar_sbst))) 
#Eliminate all values =< x, x must be varied manually
sar_sbst.r <- sar_sbst.r[,colSums(sar_sbst.r) > remove.n]
#Eliminate bottom x% of values, .x must be varied manually
sar_sbst.r <- sar_sbst.r[,rank(colSums(sar_sbst.r), ties.method = "random") > .25*(max(rank(colSums(sar_sbst.r), ties.method = "random")))]

product <-  matrix(nrow=(ncol(sar_sbst.r)),ncol=length(species),dimnames=list(colnames(sar_sbst.r),species))

for (i in species[1:length(species)]){
  product[,i] <- 1*(colSums(subset(sar_sbst.r, (mdat_sbst$Host_Spp%in%i)))>1)
  }

if(type=="euler" & length(species)>5){
product <- data.frame(product)
product <- euler(product,shape="ellipse")
return(product)
}
else if(type=="euler" & length(species)<6){
  product <- data.frame(product)
  product <- euler(product)
  return(product)
}
else if(type=="venn"){
  product <- data.frame(product)
  product <- venn(product)
  return(product)
}
}

a <- sar_venn(c("flava","rosea","leucophylla","leuco-rosea","flava-leuco","flava-rosea"))
plot(a,labels=c("Flava","Rosea","Leucophylla", "Leucophylla-Rosea", "Flava-Leucophylla","Flava-Rosea"),
           quantities=T, main="Title")

test <- sar_venn(c("rosea","leucophylla","leuco-rosea"))

#Print as Pdf
pdf("d18s_sfslsro_Venn.pdf")
plot(a, quantities=T)
dev.off()


