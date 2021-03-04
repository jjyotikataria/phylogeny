#setting up libraries
#Biostrings Bioconductor package
#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")

library(seqinr)
library(Biostrings)
library(ape)
library(textmineR)

#In the folder containing only the 12 viral sequences
nfiles <- length(dir())
seqdat <- vector("list", nfiles)

for(i in 1:nfiles){
seqdat[[i]] <- read.fasta(file=dir()[i])
}

label <- sapply(1:nfiles, function(k) unlist(strsplit(dir()[k], "[.]"))[1])
names(seqdat) <- label

seqdat_join <- lapply(seqdat, function(k) paste(toupper(unlist(k)),collapse="") )
      
The next step is to generate tabulate the frequency of all possible 5-mers in each of the viral genomes.

#enumerate all 5-mers
dna <- c("A","C","G","T")
kmer5 <- expand.grid(dna, dna, dna, dna, dna)
kmer5 <- apply(kmer5, 1, function(k) paste(k, collapse=""))

#function for counting all possible kmers (k=5) given a single dna string
kmercount <- function(data){
  sapply(1:length(kmer5), function(k)
  length(unlist(gregexpr2(kmer5[k], data)))
  )}

#vector of counts for all possible kmers (k=5) for all viral sequences
kmer_features <- lapply(seqdat_join, function(k) kmercount(k))


#Collect k-mer counts into a data frame
M <- do.call(rbind, kmer_features)
      
We then take care of some labelling stuff.

#taxonomic labels
taxonomy <- data.frame(rownames(M), 
c("Filoviridae", "Flaviviridae", "Picornaviridae", "Flaviviridae",
"Picornaviridae", "Flaviviridae", "Filoviridae", "Flaviviridae",
"Flaviviridae", "Picornaviridae", "Picornaviridae",
"Picornaviridae"))
colnames(taxonomy) <- c("Virus","Family")

#Simplify virus species names
virusnames <- sapply(1:nrow(taxonomy), function(k){
chop <- unlist(strsplit(as.character(taxonomy$Virus)[k],"_")) 
chop[length(chop)]}
)

rownames(M) <- virusnames

tipcolor <- c("red","blue","darkviolet")[unclass(taxonomy$Family)]
      
We compute the appropriate distance between all pairs of 5-mer probability vectors for the viral genomes using the Jensen-Shannon divergence, and then use the resulting distance matrix as input for the bioNJ algorithm for tree construction.

#The correct input for CalcJSDivergence is the (unnormalised) count vector
JSdist <- CalcJSDivergence(M)

plot.phylo(bionj(JSdist), type="unrooted", cex=0.8, tip.color=tipcolor,
rotate.tree=95)
      
A principal component analysis (PCA) plot is also useful for corroborating the topology of the inferred tree. Here, it seems that the plot of the first and the fourth principal component gives the best ordination, whereby viruses from the same family are clustered together.

#Alternative visualisation
#Q-mode PCA for visualisation of ordination of virus species in low dimensional kmer space
M_norm <- t(apply(M, 1, function(k) k/sum(k)))
pca <- prcomp(M_norm)
summary(pca)

pairs(pca$x[,1:5], pch=16, cex=2,
col=c("red","blue","olivedrab")[unclass(taxonomy$Family)])

#PC1 vs PC4 seems to give best separation for the three viral families

plot(pca$x[,1], pca$x[,4], cex=2, xlim=c(-0.022,0.022), ylim=c(-0.01, 0.01),
pch=16, col=c("red","blue","darkviolet")[unclass(taxonomy$Family)],
xlab="PC1(38%)", ylab="PC5(6%)")

text(pca$x[,1], pca$x[,4]+0.001, virusnames,  cex=0.7, srt=35)

legend(0.01,-0.007,pch=16, pt.cex=2, 
col=c("red","blue","darkviolet"), c("Filoviridae","Flaviviridae","Picornaviridae"))
      
