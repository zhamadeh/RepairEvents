################################################
            # Packages #
################################################
library(tidyverse)
library(GenomicRanges)

################################################
            # Data #
################################################
breakpoints <- read.table("breakPointSummary.txt",header = T)
breakpoints$filenames=tools::file_path_sans_ext(breakpoints$filenames)
breakpoints$filenames<-as.factor(breakpoints$filenames)
breakpoints$gene = NA

for (row in 1:nrow(breakpoints)){
  ID=strsplit(as.character(breakpoints[row,]$filenames),split = "[-_.]")[[1]]
  if ("blm" %in% tolower(ID) ){
    if ("recq5" %in% tolower(ID) ||"recql5" %in% tolower(ID) ){
      id <- "BLM/RECQL5"
    } else {
      id <- "BLM"
    }
  } else if ("recq5" %in% tolower(ID) ||"recql5" %in% tolower(ID) ){
    if (! "blm" %in% tolower(ID) ){
      id <- "RECQL5"
    }
  } else {id <- "WT" }
  
  breakpoints[row,]$gene=id
}

breakpoints$filenames<-as.factor(breakpoints$filenames)
breakpoints$gene<-as.factor(breakpoints$gene)


################################################
            # Summary #
################################################

summary <- breakpoints %>% group_by(filenames,gene)%>%summarize(n())
summary <- summary %>% group_by(gene)%>%summarize(n())

summaryBreaks.df=GRanges(breakpoints)
filterFrequency=

centromeres <- read.table("../StrandSeeker/Input/Blacklist/centromeres2.txt",header=F) #%>% select(-c(V4))
centromeres <-centromeres %>% dplyr::rename("seqnames"=V1,"start"=V2,"end"=V3)
centroGRange <- GRanges(centromeres)

suppressWarnings(breakGRanges <-summaryBreaks.df[-queryHits(findOverlaps(summaryBreaks.df, centroGRange, type="any")),])
breakGRanges$freq = 0

for (row in 1:length(breakGRanges)){
  tmp = breakGRanges[row,]
  query= breakGRanges[breakGRanges$gene==tmp$gene,]
  breakGRanges[row,]$freq = countOverlaps(tmp,query, maxgap = 500000)
}

breaks = merge(as.data.frame(breakGRanges),summary,by="gene")
breaks$alleleFreq=round(((breaks$freq/breaks$`n()`)*100),digits = 0)


inversions <- breakGRanges[breakGRanges$allele_freq>filterFrequency,]
breaks <- breakGRanges[breakGRanges$allele_freq<filterFrequency,]
breaks <- as.data.frame(breaks)
breaks <- filter(breaks,width<1000000)

write.table(test,"test.txt",sep = "\t",col.names = T,row.names = F,quote=F)
write.table(breaks,"breaks.txt",sep = "\t",col.names = T,row.names = F,quote=F)
recurring=filter(breaks,alleleFreq > 25)
recurring$alleleFreq<-as.factor(recurring$alleleFreq)

i=levels(recurring$gene)[1]
j=levels(tmp$alleleFreq)[1]
for (i in levels(breaks$gene)){
  tmp <- filter(recurring,gene==i)
  for (j in levels(tmp$alleleFreq)){
    tmp2 <- filter(tmp, alleleFreq==j)
    vec = tmp2$seqnames
    chr <- names(table(vec))[as.vector(table(vec))==max(table(vec))]
    tmp2 <- filter(tmp2,seqnames==chr)
    
  }
}





test <-  filter(breaks,gene=="BLM/RECQL5" & alleleFreq >57 & alleleFreq<75)


