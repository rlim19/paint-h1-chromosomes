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
unsorted_domains <- read.table("data/PyStatesWithRepeats.bed.gz")
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
  #col=hmcols<-colorRampPalette(c("green","yellow","red"))(256),
  heatmap.2(mat_, col= colorRampPalette(c("white","green","green4","violet","purple"))(100),     
            Colv=Colv, Rowv=Rowv,
            dendrogram="none", trace="none",
            density.info = 'none',scale="none",
            keysize=0.8, margins = margins,
            symbreaks=cutZero)
  
  
}

trans_domain <- getTransitionMatrix(numericDomains)

#mac heatmap
heatmap.2(trans_domain,col= colorRampPalette(c("white","yellow","red"))(80), dendrogram="none", trace="none", density.info="none")

colnames(trans_domain) <- rownames(trans_domain) <-  c('repeat','black', 'purple', 'yellow', 'red')
get_heatMat(trans_domain, margins=c(8,8), Colv=NA, Rowv=TRUE, cutZero=FALSE)
round(trans_domain*100,1)
       repeat black purple yellow pink  red
repeat    0.0  61.3   22.4   10.0  4.6  1.6
black    87.6   0.0    6.9    1.6  2.6  1.3
purple   53.1  11.4    0.0    5.7 24.2  5.6
yellow   46.3   5.0   11.3    0.0 26.8 10.5
pink     14.1   5.9   33.2   19.1  0.0 27.7
red      10.0   5.6   15.0   15.5 53.9  0.0

# calculte size's summary of each domain
# domains plus mapped with repeat domains
colorRepeat_domains <- read.table("data/PyStatesWithRepeats.bed.gz")
head(colorRepeat_domains)

domainRepeat_sizeSummary <- tapply(INDEX=colorRepeat_domains$V4, 
                      X=colorRepeat_domains$V3-colorRepeat_domains$V2,
                      summary)
domainRepeat_sizeSummary

colorRepeat_domains$size <- colorRepeat_domains$V3-colorRepeat_domains$V2
colorRepeat_size <- colorRepeat_domains[,c(4,5)]
head(colorRepeat_size)
colnames(colorRepeat_size) <- c("state", "size")
colorRepeat_size$state <- factor(colorRepeat_size$state, 
                          levels= c("black", "purple", "yellow", 
                                    "pink", "red", "repeat"))
library(ggplot2)

DomainsizePlot <- ggplot(colorRepeat_size, aes(state, size, fill=state)) + geom_boxplot(outlier.shape = NA, fatten=8) + coord_cartesian(ylim = c(0,25000)) +
scale_fill_manual(name = "Color States", 
                                   values=c("black","purple4",
                                            "gold2","deeppink2", 
                                            "red","gray" )) 

DomainsizePlot
# annotate's mean
DomainsizePlot <- DomainsizePlot+ 
#mean black's domain 15520
annotate(geom="text", x=1, y=15520, label="--------", colour="white", size=14, fontface="bold", angle=0) +
#mean purple's domain 16390
annotate(geom="text", x=2, y=16390, label="--------", colour="white", size=14, fontface="bold", angle=0) +
#mean yellow's domain 18170 
annotate(geom="text", x=3, y=18170, label="--------", colour="white", size=14, fontface="bold", angle=0) + 
#mean pink's domain 6001
annotate(geom="text", x=4, y=6001, label="--------", colour="white", size=14, fontface="bold", angle=0) +
#mean red's domain 4462
annotate(geom="text", x=5, y=4462, label="--------", colour="white", size=14, fontface="bold", angle=0) +
#mean repeat's domain 8616 
annotate(geom="text", x=6, y=8616 , label="--------", colour="white", size=14, fontface="bold", angle=0)

# enlarge the text's size
DomainsizePlot + theme(text = element_text(size=30),                                   
      axis.text.y = element_text(angle=0, vjust=1, colour="black")) + 
  opts(legend.position = "none")                                        






