################################################
            # Packages #
################################################
library(tidyverse)
library(GenomicRanges)
library(breakpointR)
library(rtracklayer)
library(doParallel)

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
breaks$alleleFreq=round((breaks$freq/breaks$n..*100),digits = 0)

breaks$gene=str_replace(breaks$gene,"/","-")
breaks$gene <- as.factor(breaks$gene)

recurring=filter(breaks,alleleFreq > 25)
recurring$alleleFreq<-as.factor(recurring$alleleFreq)

recurring$gene <- as.factor(recurring$gene)
recurring$alleleFreq<-as.factor(recurring$alleleFreq)

write.table(breaks,"DATA/breaks.txt",quote = F,row.names = F,col.names = T,sep="\t")
write.table(recurring,"DATA/recurringBreaksCutoff25.txt",quote = F,row.names = F,col.names = T,sep="\t")

breaks <- read.table("DATA/breaks.txt",header=T)
breaks <- read.table("DATA/breaks.txt",header=T)


################################################
                  # functions #
################################################


repairER <- function(i, recurring, breaks){
  tmp <- filter(recurring,gene==i)
  for (j in levels(tmp$alleleFreq)){
    tmp2 <- filter(tmp, alleleFreq==j)
    tmp2$seqnames<-droplevels(tmp2$seqnames)
    for (t in levels(tmp2$seqnames)){
      chr=t
      tmp3 <- filter(tmp2,seqnames==chr)
      if (nrow(tmp3) > 5){
        message("Found a possible inversion on ", chr," in ", j,"% of ",i," cells.")
        
        tmp3$filenames <- droplevels(tmp3$filenames)
        
        plotspath=paste0("PLOTS/",chr,"-",i,"-",j,"-f/")
        datapath=paste0(plotspath,"data/")
        
        if (!file.exists(plotspath) ) { dir.create(plotspath)}
        if (!file.exists(datapath) ) { dir.create(datapath)}
        
        #use tmp2 so you get other breakpoints that were also on the same chromosome
        tmp4 <- filter(breaks,gene==i)
        bed <- tmp4 %>% filter(seqnames==chr & filenames %in% tmp3$filenames) %>% select(seqnames,start,end,filenames)
        bedpath=paste0(plotspath,"breakpoints.bed")
        export(bed,bedpath,format = "gff3")
        
        readsFiles<- tools::file_path_sans_ext(levels(tmp3$filenames))
        files2transfer=paste0("DATA/browserfiles/",readsFiles,"_reads.bed.gz")
        
        file.copy(files2transfer,datapath)
        
        files2plot=paste0("DATA/browserfiles/",levels(tmp3$filenames))
        plotBreakpointsPerChr(files2plot,plotspath = plotspath,chromosomes = c(chr))
      }
    }
  }
}

startTimedMessage <- function(...) {
  
  x <- paste0(..., collapse='')
  message(x, appendLF=FALSE)
  ptm <- proc.time()
  return(ptm)
  
}


stopTimedMessage <- function(ptm) {
  
  time <- proc.time() - ptm
  message(" ", round(time[3],2), "s")
  
}

################################################
     # Running functions generate output #
################################################

ptm <- startTimedMessage("Reading file ", bamfile, " ...")

for (i in levels(recurring$gene)){
  tmp <- filter(recurring,gene==i)
  for (j in levels(tmp$alleleFreq)){
    tmp2 <- filter(tmp, alleleFreq==j)
    tmp2$seqnames<-droplevels(tmp2$seqnames)
    for (t in levels(tmp2$seqnames)){
      chr=t
      tmp3 <- filter(tmp2,seqnames==chr)
      if (nrow(tmp3) > 5){
        message("Found a possible inversion on ", chr," in ", j,"% of ",i," cells.")
        
        tmp3$filenames <- droplevels(tmp3$filenames)
        
        plotspath=paste0("PLOTS/",chr,"-",i,"-",j,"-f/")
        datapath=paste0(plotspath,"data/")
        
        if (!file.exists(plotspath) ) { dir.create(plotspath)}
        if (!file.exists(datapath) ) { dir.create(datapath)}
        
        #use tmp2 so you get other breakpoints that were also on the same chromosome
        tmp4 <- filter(breaks,gene==i)
        bed <- tmp4 %>% filter(seqnames==chr & filenames %in% tmp3$filenames) %>% select(seqnames,start,end,filenames)
        bedpath=paste0(plotspath,"breakpoints.bed")
        export(bed,bedpath,format = "gff3")
        
        readsFiles<- tools::file_path_sans_ext(levels(tmp3$filenames))
        files2transfer=paste0("DATA/browserfiles/",readsFiles,"_reads.bed.gz")
        
        file.copy(files2transfer,datapath)
        
        files2plot=paste0("DATA/rSdata/",levels(tmp3$filenames))
        plotBreakpointsPerChr(files2plot,plotspath = plotspath,chromosomes = c(chr))
      }
    }
  }
}

stopTimedMessage(ptm)



ptm <- startTimedMessage()

coresToUse = detectCores()/2
cl =makeCluster(coresToUse)
doParallel::registerDoParallel(cl)
foreach::foreach(i=levels(recurring$gene),.packages = c("rtracklayer","tidyverse","GenomicRanges","breakpointR"),.export = c("repairER")) %dopar% {repairER(i,recurring,breaks)}
stopCluster(cl)

stopTimedMessage(ptm)



