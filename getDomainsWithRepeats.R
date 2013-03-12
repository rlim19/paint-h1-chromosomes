###############################
# get the coordinates of NAs  #
###############################

# All_bigTable containing all the NAs and 0s from 93 Profiles mapped with GEM
data <- read.delim("data/ColorAll_bigTable.txt")
head(data)
data$tag <- as.character(data$tag)
dim(data)

# get the rows for repeats
row_na1 <- which(rowSums(is.na(data[,5:97])) > 0)
row_na2 <- which(rowSums(data[,5:97]) == 0)
row_na <- c(row_na1, row_na2)
head(row_na1)
head(row_na2)
length(row_na)
head(row_na)

data[row_na, "tag"] <- "repeat"
head(data)

# create bed format
data_bed <- data[,c(2,3,4,1)]
dim(data_bed)
head(data_bed)

# write the domains (including repeats) in bed 
write.table(data_bed, "data/ExpandedStatesWithRepeats_dom_hESC.bed", quote=FALSE, row.names=FALSE, col.names=FALSE,
            sep="\t")

getDomain <- function(data_f, domain_name){
  
  domain <- data_f
  domain$tag <- domain$tag==as.character(domain_name)
  domain <- domainify(domain)
  domain$tag <- rep(domain_name, nrow(domain))
  return(domain)
  
}

# compressed the bed format into domains(concatenating consecutive states)
source("domainify.R")
All_domains = data.frame()
for (i in levels(as.factor(data_bed$tag))){
 dom <- getDomain(data_f=data_bed, domain_name=i)
 All_domains <- rbind(All_domains, dom)
}

write.table(All_domains, "data/statesWithRepeats_dom_hESC.bed", quote=FALSE, row.names=FALSE, col.names=FALSE,
            sep="\t")

# check the transition matrix of domains
unsorted_domains <- read.table("data/PyStatesWithRepeats.bed")
head(unsorted_domains)
dim(unsorted_domains)
source('/users/gfilion/rlim/R_misc/getTransitionMatrix.R')
c('black', 'purple', 'yellow', 'pink','red', 'repeat'))

numericDomains <- match(unsorted_domains$V4, c('repeat','black', 'purple', 'yellow', 'pink', 'red'))
length(numericDomains)

library(gplots)
get_heatMat <- function(mat_, margins, Colv, Rowv, cutZero){
  # construct heatmap, given a matrix
  
  library(gplots)
  heatmap.2(mat_, col=hmcols<-colorRampPalette(c("green","yellow","red"))(256), 
            Colv=Colv, Rowv=Rowv,
            dendrogram="none", trace="none",
            density.info = 'none',scale="none",
            keysize=0.8, margins = margins,
            symbreaks=cutZero)
  
  
}

trans_domain <- getTransitionMatrix(numericDomains)


colnames(trans_domain) <- rownames(trans_domain) <-  c('repeat','black', 'purple', 'yellow', 'pink', 'red')
get_heatMat(trans_domain, margins=c(5,10), Colv=NA, Rowv=TRUE, cutZero=FALSE)
round(trans_domain*100,1)
       repeat black purple yellow pink  red
repeat    0.0  61.3   22.4   10.0  4.6  1.6
black    87.6   0.0    6.9    1.6  2.6  1.3
purple   53.1  11.4    0.0    5.7 24.2  5.6
yellow   46.3   5.0   11.3    0.0 26.8 10.5
pink     14.1   5.9   33.2   19.1  0.0 27.7
red      10.0   5.6   15.0   15.5 53.9  0.0