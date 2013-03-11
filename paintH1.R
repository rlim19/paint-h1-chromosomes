####################################################### 
# paint chromosomes in H1 with color chromatin states #
#######################################################

library(ggbio)
library(rtracklayer)
# import color states data
h1_state <- import("data/states_dom_hESC.bed.gz")
head(h1_state)
colnames(h1_state) <- "State"
h1_state <- as(h1_state, "GRanges")

library(BSgenome.Hsapiens.UCSC.hg19)
chr.len = seqlengths(Hsapiens)
chr.len
#exclude chromosomes with suffix "_" , "M", "Het", "extra".
chr.len = chr.len[grep("_|M|U|Het|extra", names(chr.len), invert = T)] 
#order the chromosomes
h1_state= keepSeqlevels(h1_state, names(chr.len))   
seqlevels(h1_state) = names(chr.len) 
seqlengths(h1_state) = (chr.len)
print(h1_state)

chrom.col <- c("black", "purple4", "gold2","deeppink2", "red")
h1_state$State <- factor(h1_state$State, 
                         levels=c('black', 'purple', 'yellow', 'pink','red'))
p <- autoplot(h1_state, layout = "karyogram", aes(fill = State))
p + scale_fill_manual(values = chrom.col) + opts(legend.position = "none") 


# circos's style
p <- autoplot(h1_state, layout = "circle", aes(fill = State))
