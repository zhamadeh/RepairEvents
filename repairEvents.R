################################################
            # Packages #
################################################
library(tidyverse)
library(GenomicRanges)
library(breakpointR)

################################################
            # Data #
################################################
breakpoints <- read.table("breakPointSummary.txt",header = T)
breakpoints$filenames=tools::file_path_sans_ext(breakpoints$filenames)
breakpoints$filenames<-as.factor(breakpoints$filenames)
breakpoints$gene = NA

#ANNOTATE with gene information (slow)
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

#Set variables as factors
breakpoints$filenames<-as.factor(breakpoints$filenames)
breakpoints$gene<-as.factor(breakpoints$gene)


################################################
            # Summary #
################################################

sum <- breaks %>% group_by(filenames,gene)%>%summarize(n())
summary <- sum %>% group_by(gene)%>%summarize(n())

summaryBreaks.df=GRanges(breakpoints)

#Blacklist centromeress
centromeres <- read.table("../StrandSeeker/Input/Blacklist/centromeres2.txt",header=F) #%>% select(-c(V4))
centromeres <-centromeres %>% dplyr::rename("seqnames"=V1,"start"=V2,"end"=V3)
centroGRange <- GRanges(centromeres)

suppressWarnings(breakGRanges <-summaryBreaks.df[-queryHits(findOverlaps(summaryBreaks.df, centroGRange, type="any")),])


################################################
                # Allele Frequency #
################################################

breakGRanges$freq = 0

for (row in 1:length(breakGRanges)){
  tmp = breakGRanges[row,]
  query= breakGRanges[breakGRanges$gene==tmp$gene,]
  breakGRanges[row,]$freq = countOverlaps(tmp,query, maxgap = 500000)
}

breaks = merge(as.data.frame(breaks),summary,by="gene")
breaks$alleleFreq=round(((breaks$freq/breaks$`n()`)*100),digits = 0)


#breaks <- read.table("breaks.txt",header=T)





################################################
                  # plotting #
################################################

inversions <- breakGRanges[breakGRanges$allele_freq>filterFrequency,]
breaks <- breakGRanges[breakGRanges$allele_freq<filterFrequency,]
breaks <- as.data.frame(breaks)

breaks <- filter(breaks,width<1000000)

breaks$gene=str_replace(breaks$gene,"/","-")
breaks$gene <- as.factor(breaks$gene)
breaks$alleleFreq<-as.factor(breaks$alleleFreq)

recurring=filter(breaks,alleleFreq > 25)
recurring$alleleFreq<-as.factor(recurring$alleleFreq)

recurring$gene <- as.factor(recurring$gene)
recurring$alleleFreq<-as.factor(recurring$alleleFreq)

for (i in levels(recurring$gene)){
  tmp <- filter(recurring,gene==i)
  for (j in levels(tmp$alleleFreq)){
    tmp2 <- filter(tmp, alleleFreq==j)
    tmp2$seqnames<-droplevels(tmp2$seqnames)
    for (t in levels(tmp2$seqnames)){
      chr=t
      tmp3 <- filter(tmp2,seqnames==chr)
      if (nrow(tmp3) > 5){
        tmp3$filenames <- droplevels(tmp3$filenames)
        
        bed <- tmp3 %>% select(seqnames,start,end,filenames)
        export(bed,"bed.bed",format = "gff3")
        
        plotspath=paste0("PLOTS/",chr,"-",i,"-",j,"-f/")
        datapath=paste0(plotspath,"data/")
        
        if (!file.exists(plotspath) ) { dir.create(plotspath)}
        if (!file.exists(datapath) ) { dir.create(datapath)}
        
        
        readsFiles<- tools::file_path_sans_ext(levels(tmp3$filenames))
        files2transfer=paste0("../data/browserfiles/",readsFiles,"_reads.bed.gz")
        
        file.copy(files2transfer,datapath)
        
        files2plot=paste0("../data/data/",levels(tmp3$filenames))
        plotBreakpointsPerChr(files2plot,plotspath = plotspath,chromosomes = c(chr))
      }
    }
  }
}


