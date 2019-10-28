###2015 Sarracenia Data Venn Diagram Comparisons###

##Load relevant packages
if (!require("eulerr")) {install.packages("eulerr"); require("eulerr")}
if (!require("vegan")) {install.packages("vegan"); require("vegan")}

##Load Data
setwd("workingdirectory_path")
#Raw Data
sar <- read.table("github_tutorial_data.txt", header = T, row.names = 1, check.names = F, sep = "\t") 
sar <- sar[order(row.names(sar)),]
#Metadata
mdat <- read.csv("github_tutorial_metadata.csv", head=T, row.names=1) 
mdat <- mdat[order(row.names(mdat)),]
#Remove shotgun samples
sar <- subset(sar, mdat$Project %in% ("Comparison"))
mdat <- subset(mdat, row.names(mdat) %in% row.names(sar))
#Check that tables match
row.names(mdat) == row.names(sar)

###Function: sar_venn()
##Syntax: sar_venn(species,type,remove.n)
##Arguments: 
#species = any combination of "alata", "flava", "leucophylla", "purpurea", "rosea", "rubra", "flava-leuco", "flava-rosea", "leuco-rosea"
#type = character, "euler" or "venn", default="euler"
#remove.n = numeric, default=0, removes x number of OTUs off the bottom

###Start function
sar_venn <- function(species,type="euler",remove.n=0){
  
#Subset raw data of given species
sar_sbst <- subset(sar, mdat$Host_Spp %in% species)
#Subset metadata of given species
mdat_sbst <- subset(mdat, row.names(mdat) %in% row.names(sar_sbst))
#Drop old factors from metadata
mdat_sbst$Host_Spp <- mdat_sbst$Host_Spp[drop=T]

#Rarify
sar_sbst.r <- rrarefy(sar_sbst, min(rowSums(sar_sbst))) 
#Eliminate all values =< remove.n
sar_sbst.r <- sar_sbst.r[,colSums(sar_sbst.r) > remove.n]

#create empty dataframe to store presence/absence vectors
product <-  matrix(nrow=(ncol(sar_sbst.r)),ncol=length(species),dimnames=list(colnames(sar_sbst.r),species))

#create vector for presence/absence of OTUs for each given species
for (i in species[1:length(species)]){ 
  product[,i] <- 1*(colSums(subset(sar_sbst.r, (mdat_sbst$Host_Spp%in%i)))>1)
  }
  
#generate diagrams
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

#test
a <- sar_venn(c("rosea","leucophylla","leuco-rosea"))
plot(a, quantities=T)

#Print as Pdf
pdf("test_venn.pdf")
plot(a, quantities=T)
dev.off()


