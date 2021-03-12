################################################
              # Packages #
################################################
library(tidyverse)
library(GenomicRanges)
library(breakpointR)
library(rtracklayer)
library(doParallel)

################################################
            # Core Functions #
################################################

#Collect all breakpoints from RData files in datapath

collectBreaksAllFiles <- function(datapath="DATA/rdata/"){
  files <- list.files(datapath, pattern=".RData$", full.names=TRUE)
  
 
  breaks.all.files <- list()
  breaksConfInt.all.files <- list()
  summaryBreaks <- list()
  n=1
  for (file in files) {
    message("Reading ... " , basename(file), " ... ",round(  (n/length(files))*100  ,  digits = 1  ) , "%"  )
    n=n+1
    data <- get(load(file))[c('breaks', 'confint','ID')]
    data$breaks$ID <- data$ID
    summaryBreaks[[basename(file)]] <- summarizeBreaks(data)
    breakpoints <- data$breaks
    breaks.confint <- data$confint
    if (length(breakpoints)) {
      suppressWarnings( breaks.all.files[[file]] <- breakpoints ) #TODO check if this can be done without warnings
    }  
    if (length(breaks.confint)) {
      suppressWarnings( breaksConfInt.all.files[[file]] <- breaks.confint ) 
    }  
  }
  return(breaks.all.files)
}


#Use p-values to call hotspots for breakpoints

breakpointHotspotter <- function(breaks.all.files){
  
  #Take in all breakpoints and calcualte p-values
  
  names(breaks.all.files) <- NULL
  gr.list=breaks.all.files
  bw=1000000
  pval=1e-8
  names(gr.list) <- NULL
  gr <- do.call(c, gr.list)
  gr <- GenomicRanges::sort(gr)
  
  ## Iterate over chromosomes and calculate p-values
  pranges.list <- GenomicRanges::GRangesList()
  hotspots=data.frame()
  count=1
  
  df = data.frame(chr=c(),midpoint=c(),density=c(),pvalue=c(),null_midpoints=c(),null_density=c())
  df_p=data.frame(pval=c(),chrom=c())
  
  for (chrom in seqlevels(gr)) {
    grc <- gr[seqnames(gr)==chrom]
    if (length(grc)>1) {
      midpoints <- (start(grc)+end(grc))/2
      kde <- stats::density(midpoints,bw=bw,kernel='gaussian')
      # Random distribution of genomic events
      kde.densities <- numeric()
      
      for (i1 in seq_len(100)) {
        midpoints.r <- round(stats::runif(length(midpoints),1,seqlengths(gr)[chrom]))
        kde.r <- stats::density(midpoints.r,bw=bw,kernel='gaussian')
        kde.densities <- c(kde.densities, kde.r$y)
      }
      # Use ecdf to calculate p-values 
      p <- 1-stats::ecdf(kde.densities)(kde$y)
      pvalues <- data.frame(chromosome=chrom,start=kde$x,pvalue=p)
      # Make GRanges
      pvalues$end <- pvalues$start
      pvalues$chromosome <- factor(pvalues$chromosome, levels=seqlevels(gr))
      pvalues <- as(pvalues,'GRanges')
      seqlevels(pvalues) <- seqlevels(gr)
      suppressWarnings(
        seqlengths(pvalues) <- seqlengths(gr)[names(seqlengths(pvalues))]
      )
      # Resize from pointsize to bandwidth
      suppressWarnings(
        pvalues <- GenomicRanges::resize(pvalues, width=bw, fix='center')
      )
      pvalues <- trim(pvalues)
      ## Find regions where p-value is below specification
      mask <- pvalues$pvalue <= pval
      rle.pvals <- rle(mask)
      rle.pvals$values <- cumsum(rle.pvals$values+1)
      pvalues$group <- inverse.rle(rle.pvals)
      if (length(which(mask))>0) {
        pvalues.split <- split(pvalues[mask],pvalues$group[mask])
        pranges <- unlist(endoapply(pvalues.split, function(x) { y <- x[1]; end(y) <- end(x)[length(x)]; y$pvalue <- min(x$pvalue); return(y) }))
        pranges$group <- NULL
        pranges$num.events <- GenomicRanges::countOverlaps(pranges,grc)
        pranges.list[[chrom]] <- pranges
      }
      for (el in 1:length(pranges)){
        tmp = as.data.frame(gr[queryHits(findOverlaps(gr,pranges[el],type = "any"))])
        tmp$count=count
        hotspots <- rbind(tmp,hotspots)
        count=count+1
      }
      
      
      
    mp = kde$x
    ds = kde$y
    pv = pvalues$pvalue
    exp = data.frame(mp=mp,ds=ds)
    exp$type="BREAKPOINTS"
    exp$p = pv
    
    null_mp = kde.r$x
    null_ds = kde.r$y
    null = data.frame(mp=null_mp,ds=null_ds)
    null$type = "NULL"
    null$p = pv
    
    tmp = rbind(exp,null)
    tmp$chr = chrom

    #tmp2 = data.frame(pvalues=pv)
    #tmp2$chr = chrom
    
    df <- rbind(tmp,df)
    df_p <- rbind(tmp2,df_p)
    }
    count=count+1
  }
  pranges <- unlist(pranges.list, use.names=FALSE)
  names(pranges) <- NULL
  

 return(hotspots)
}


# Create file structure for each inversion: 
  ## reads.bed files for all libraries involved
  ## breakpoints.bed for breakpoint coordinates
  ## chr_breakpoints.pdf for chromosome-specific ideogram plotting
  ## genotype.txt for summary of genotype results

# One summary file for all inversions, each row:
  ## chr, mean(start), mean(end), mean(width), # of libraries, % BLM, % RECQL5, % BLM-RECQL5, % WT 



savingAndPrinting <- function(hotspots,hotpath="HOTSPOT_EVENTS",printing=F,export=F,normalize=F){
  
  # Directory for creating file structure
  
  if (!file.exists(hotpath) ) { dir.create(hotpath)}
  
  # Count corresponds to a unique inversion so convert to factor for easier iterating
  hotspots$count=as.factor(hotspots$count)
  # Also easier for dealing with IDs as factor
  hotspots$ID <- as.factor(hotspots$ID)
  
  numOfLibsPerGene<-read.table("DATA/numOfLibsPerGene.txt",header=T)
  
  # Initialize empty dataframe for summary
  summary <- data.frame(chr=c(),start= c(),end=c(),count=c(), width=c(),n=c(),perc=c(),BLM=c(),RECQL5=c(),BLM_RECQL5=c(),WT=c())
  
  # Iterate through each level of count (unique inversion) and 1) filter 2) summarize 3) save 4) plot
  hotspots$count<-as.factor(hotspots$count)
  
  for (level in levels(hotspots$count)){
    
    # 1) FILTER
    tmp <- filter(hotspots, count==level)
    chr = as.character(tmp[1,1]) # Set chromosome
    
    tmp$ID <- droplevels(tmp$ID) # Drop unused levles
    tmp$gene = NA
    for (i in 1:length(tmp$ID)){
      ID <- strsplit(as.character(tmp$ID[i]),split = "[-_]")[[1]]
      for (j in ID){
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
      }
      tmp[i,]$gene=id
    }
    perc <- tmp %>% group_by(gene) %>% summarize(perc = round((n()/nrow(tmp))*100,digits = 1),n=n())
    perc$resolution = mean(tmp$width)
    
    
    j=round((nrow(tmp)/317)*100,2)
    message("Found hotspot #",level," on ", chr," in ", nrow(tmp)," cells (",j,"%)")
    
    
    if (perc[which.max(perc$perc),]$perc > 90){
      message("This inversion is unique to ", perc[which.max(perc$perc),]$gene, " making up ",perc[which.max(perc$perc),]$perc, " of the libraries")
    }
    
    
    datapath=paste0(hotpath,"/",chr,"-",j,"/")
    if (!file.exists(datapath) ) { dir.create(datapath)}
    readsdatapath=paste0(hotpath,"/",chr,"-",j,"/reads/")
    if (!file.exists(readsdatapath) ) { dir.create(readsdatapath)}
    
    if (export){
      bed=tmp
      export(bed,paste0(datapath,"breakpoints.bed"),format = "gff3")
      write.table(as.data.frame(perc),paste0(datapath,"genotype.txt"),row.names = F,col.names = T,quote = F,sep="\t")
      export(bed,paste0(datapath,"breakpoints.bed"),format = "bed")
    }
    
    if (length(tmp$ID)>25){
      files2transfer=paste0("DATA/browserfiles/",tmp$ID[1:25],"_reads.bed.gz")
    } else {
      files2transfer=paste0("DATA/browserfiles/",tmp$ID,"_reads.bed.gz")
    }

    
    numLibs=tmp %>% group_by(gene)%>%summarize(numLibs=length(levels(droplevels(ID))))
    if(normalize==FALSE){
      if ("BLM" %in% numLibs$gene){
        BLM=as.numeric((numLibs[numLibs$gene=="BLM",2]))
      } else { BLM=0 }
      if ("RECQL5" %in% numLibs$gene){
        RECQL5=as.numeric((numLibs[numLibs$gene=="RECQL5",2]))
      } else { RECQL5=0 }
      if ("BLM/RECQL5" %in% numLibs$gene){
        BLM_RECQL5=as.numeric((numLibs[numLibs$gene=="BLM/RECQL5",2]))
      } else { BLM_RECQL5=0 }
      if ("WT" %in% numLibs$gene){
        WT=as.numeric((numLibs[numLibs$gene=="WT",2]))
      } else { WT=0 }
    } else if(normalize=="By_Library"){
      if ("BLM" %in% numLibs$gene){
        BLM=as.numeric((numLibs[numLibs$gene=="BLM",2])/numOfLibsPerGene[numOfLibsPerGene$gene=="BLM",2])
      } else { BLM=0 }
      if ("RECQL5" %in% numLibs$gene){
        RECQL5=as.numeric((numLibs[numLibs$gene=="RECQL5",2])/numOfLibsPerGene[numOfLibsPerGene$gene=="RECQL5",2])
      } else { RECQL5=0 }
      if ("BLM/RECQL5" %in% numLibs$gene){
        BLM_RECQL5=as.numeric((numLibs[numLibs$gene=="BLM/RECQL5",2])/numOfLibsPerGene[numOfLibsPerGene$gene=="BLM/RECQL5",2])
      } else { BLM_RECQL5=0 }
      if ("WT" %in% numLibs$gene){
        WT=as.numeric((numLibs[numLibs$gene=="WT",2])/numOfLibsPerGene[numOfLibsPerGene$gene=="WT",2])
      } else { WT=0 }
    } else  if (normalize=="to_100"){
      if ("BLM" %in% perc$gene){
        BLM=as.numeric((perc[perc$gene=="BLM",2]))
      } else { BLM=0 }
      if ("RECQL5" %in% perc$gene){
        RECQL5=as.numeric((perc[perc$gene=="RECQL5",2]))
      } else { RECQL5=0 }
      if ("BLM/RECQL5" %in% perc$gene){
        BLM_RECQL5=as.numeric((perc[perc$gene=="BLM/RECQL5",2]))
      } else { BLM_RECQL5=0 }
      if ("WT" %in% perc$gene){
        WT=as.numeric((perc[perc$gene=="WT",2]))
      } else { WT=0 }
    }
    
    
    row <- data.frame(chr=chr,start= mean(tmp$start),end=mean(tmp$end),count=tmp$count[1], width=mean(tmp$width),n=nrow(tmp),perc=j,BLM=BLM,RECQL5=RECQL5,BLM_RECQL5=BLM_RECQL5,WT=WT)
    summary <- rbind(summary,row)
    
    if (printing){
      file.copy(files2transfer,readsdatapath)
      files2plot=paste0("DATA/rdata/",levels(tmp$ID),".RData")
      plotBreakpointsPerChr(files2plot,plotspath = datapath,chromosomes = c(chr))
    }
  }
  return(summary)
}



################################################
        # Running Core Functions #
################################################

master <- function(){
  breaks.all.files <- collectBreaksAllFiles(datapath="DATA/rdata/")
  
  hotspots <- breakpointHotspotter(breaks.all.files)
  
  summary <- savingAndPrinting(hotspots,hotpath="HOTSPOT_EVENTS/",printing = T,normalize="By_Library",export = T)
  #summary <- savingAndPrinting(hotspots,hotpath="HOTSPOT_EVENTS/",printing = F,normalize=FALSE,export = F)
  #summary <- savingAndPrinting(hotspots,hotpath="HOTSPOT_EVENTS/",printing = F,normalize="to_100")
  
  message("\nI found ",nrow(summary), " inversions.\n")
  message("\nThe average resolution is ", round(mean(summary$width),digits = -3),".\n")
  
  write.table(summary,"summary.txt",col.names = T,row.names = F, quote = F,sep="\t")
}

summary=master()

################################################
            # Summary #
################################################

# THIS IS FOR PLOTTING 

#PREPROCESS
max_xy <-max(df$ds)
df$mp_Mb = df$mp/1e+06
df$chr <- factor(df$chr,levels = c("chr1" , "chr2",  "chr3",  "chr4" , "chr5" , "chr6",  "chr7",  "chr8" , "chr9", "chr10" ,"chr11" ,"chr12" ,"chr13", "chr14", "chr15" ,"chr16" ,"chr17", "chr18","chr19" ,"chr20", "chr21" ,"chr22", "chrX" ))

#PLOTTING
ggplot(df)+geom_point(aes(mp_Mb,ds,color=type),size=1)+facet_wrap(~chr,scales="free")+
  #geom_point(aes(mp,p))+facet_wrap(~chr,scales="free")+
  scale_y_continuous(limits = c(0,max_xy)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.title =element_text(size=28),
        text = element_text(size=20),
        legend.position = c(0.7, 0.07))+
  labs(y="DENSITY",x="CHROMOSOME POSITION")+
  guides(color = guide_legend(override.aes = list(size=10)))+
  ggsave("PLOTS/genomicDistributionOfEvents_allChr.png")

ggplot(filter(df,chr=="chr1"))+geom_point(aes(mp_Mb,ds,color=type),size=5)+facet_wrap(~chr,scales="free")+
  #geom_point(aes(mp,p))+facet_wrap(~chr,scales="free")+
  #scale_y_continuous(limits = c(0,max_xy)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.title =element_text(size=28),
        text = element_text(size=20),
        legend.position = c(0.7, 0.3))+
  labs(y="DENSITY",x="CHROMOSOME POSITION")+
  guides(color = guide_legend(override.aes = list(size=10)))+
  ggsave("PLOTS/genomicDistributionOfEvents_chr1.png")

  

ggplot(df_p)+geom_density(aes(pvalues))+facet_wrap(~chr,scales="free")






summary$count<- seq(1:nrow(summary))
summarize <- gather(summary, gene,freq,c(BLM,RECQL5,BLM_RECQL5,WT))
summarize$chr<-factor(summarize$chr,levels=c("chr1"  ,"chr2" ,"chr3" , "chr4"  ,"chr5" , "chr7",  "chr8" , "chr9"  ,"chr10", "chr11" ,"chr12" ,"chr13", "chr14"
                      ,"chr15" ,"chr16", "chr17", "chr18", "chr19" , "chr20", "chr21", "chr22", "chrX"))

summarize = summarize[order(summarize$chr),]
summarize$cat= as.factor(paste0(summarize$chr,"-",summarize$count))

summarize$cat = factor(summarize$cat,levels=c( "chr1-1" ,"chr1-2", "chr2-3" ,  "chr2-4" ,  "chr3-5" ,  "chr4-6" ,  "chr5-7" ,  "chr5-8"  , "chr7-10",  "chr7-11" , "chr7-9"   ,"chr8-12" 
                                               ,"chr8-13" , "chr9-14" , "chr9-15" , "chr9-16" 
                                              , "chr10-17" , "chr11-18", "chr11-19" ,  "chr12-20", "chr13-21", "chr14-22", "chr15-23" ,"chr15-24" ,"chr15-25"
                                             ,"chr15-26" ,"chr16-27", "chr16-28" ,"chr16-29", "chr17-30", "chr18-31", "chr19-32", "chr19-33", "chr20-34" ,"chr21-35" ,"chr22-36"
                                             ,"chr22-37" ,"chrX-38"  ,"chrX-39" ))

ggplot(summarize)+geom_col(aes(cat,freq,fill=gene))+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  labs(y="Number of events/library",x="Hotspots")
ggplot(summarize)+geom_col(aes(chr,freq,fill=gene))+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  labs(y="Number of events",x="Hotspots")
summarize$chr <- as.factor(summarize$chr)







cfs <-read.table("DATA/CFS_human.bed",fill=T) %>% select(c(V1,V2,V3))
colnames(cfs)<- c("chr","start","end")
cfs$type="CFS"
cfs$chr=factor(cfs$chr,levels=c("chr1"  ,"chr2" ,"chr3" , "chr4"  ,"chr5" ,"chr6", "chr7",  "chr8" , "chr9"  ,"chr10", "chr11" ,"chr12" ,"chr13", "chr14"
                                      ,"chr15" ,"chr16", "chr17", "chr18", "chr19" , "chr20", "chr22", "chrx"))

sum=select(summary,c(chr,start,end))
sum$type="BREAK_HOTSPOT"
sum$chr=factor(sum$chr,levels=c("chr1"  ,"chr2" ,"chr3" , "chr4"  ,"chr5" , "chr7",  "chr8" , "chr9"  ,"chr10", "chr11" ,"chr12" ,"chr13", "chr14"
                                      ,"chr15" ,"chr16", "chr17", "chr18", "chr19" , "chr20", "chr21", "chr22", "chrX"))

s=rbind(sum,cfs)
s = as.data.frame(s %>% group_by(type,chr)%>% summarize(n=n()))
s[44,2]="chrX"
s$chr=factor(s$chr,levels=c("chr1"  ,"chr2" ,"chr3" , "chr4"  ,"chr5" ,"chr6", "chr7",  "chr8" , "chr9"  ,"chr10", "chr11" ,"chr12" ,"chr13", "chr14"
                                      ,"chr15" ,"chr16", "chr17", "chr18", "chr19" , "chr20","chr21", "chr22", "chrx","chrX"))

ggplot()+ geom_bar(s,mapping=aes(chr,n,fill=type),stat="identity",position=  position_stack(reverse=T))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  labs(x="Chromosome",y="Count")





