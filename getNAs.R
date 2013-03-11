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

getDomain <- function(data_f, domain_name){
  
  domain <- data_f
  domain$tag <- domain$tag==as.character(domain_name)
  domain <- domainify(domain)
  domain$tag <- rep(domain_name, nrow(domain))
  return(domain)
  
}
source("domainify.R")
All_domains = data.frame()
for (i in levels(as.factor(data_bed$tag))){
 dom <- getDomain(data_f=data_bed, domain_name=i)
 All_domains <- rbind(All_domains, dom)
}

head(All_domains)

write.table(All_domains, "data/statesWithRepeats_dom_hESC.bed", quote=FALSE, row.names=FALSE, col.names=FALSE,
            sep="\t")
