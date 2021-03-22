rdat################################################
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
    #df_p <- rbind(tmp2,df_p)
    }
    count=count+1
  }
  pranges <- unlist(pranges.list, use.names=FALSE)
  names(pranges) <- NULL
  
  write.table(df,"SUMMARY/densityPvalueSummary.txt",col.names = T,row.names = F, quote = F,sep="\t")
  
  return(hotspots)
}


# Create file structure for each inversion: 
  ## reads.bed files for all libraries involved
  ## breakpoints.bed for breakpoint coordinates
  ## chr_breakpoints.pdf for chromosome-specific ideogram plotting
  ## genotype.txt for summary of genotype results

# One summary file for all inversions, each row:
  ## chr, mean(start), mean(end), mean(width), # of libraries, % BLM, % RECQL5, % BLM-RECQL5, % WT 



savingAndPrinting <- function(hotspots,hotpath="HOTSPOT_EVENTS",printing=F,export=F,normalize=F,cfs=F,genomeInstability=T){
  
  # Directory for creating file structure
  
  if (!file.exists(hotpath) ) { dir.create(hotpath)}
  
  # Count corresponds to a unique inversion so convert to factor for easier iterating
  hotspots$count=as.factor(hotspots$count)
  # Also easier for dealing with IDs as factor
  hotspots$ID <- as.factor(hotspots$ID)
  
  numOfLibsPerGene<-read.table("../SV_hotspot_backupCopies/Structural_Variant_Hotspotter-2/DATA/numOfLibsPerGene.txt",header=T)
  
  # Initialize empty dataframe for summary
  summary <- data.frame(chr=c(),start= c(),end=c(),count=c(), width=c(),n=c(),perc=c(),BLM=c(),RECQL5=c(),BLM_RECQL5=c(),WT=c())
  
  # Iterate through each level of count (unique inversion) and 1) filter 2) summarize 3) save 4) plot
  hotspots$count<-as.factor(hotspots$count)
  
  for (level in levels(hotspots$count)){
    message("Level: ",level)
    
    # 1) FILTER
    tmp <- filter(hotspots, count==level)
    chr = as.character(tmp[1,1]) # Set chromosome
    
    tmp$ID <- droplevels(tmp$ID) # Drop unused levles
    tmp$gene = NA
    
    if (genomeInstability==T){
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
    }
    perc <- tmp %>% group_by(gene) %>% summarize(perc = round((n()/nrow(tmp))*100,digits = 1),n=n())
    perc$resolution = mean(tmp$width)
    
    
    j=round((length(levels(droplevels(tmp$ID)))/317)*100,2)
    message("Found hotspot #",level," on ", chr," in ", nrow(tmp)," cells (",j,"%)")
    
    if (genomeInstability==T){
      if (perc[which.max(perc$perc),]$perc > 90){
        message("This inversion is unique to ", perc[which.max(perc$perc),]$gene, " making up ",perc[which.max(perc$perc),]$perc, " of the libraries")
      }
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

    if (genomeInstability==T){
      numLibs=tmp %>% group_by(gene)%>%summarize(numLibs=length(levels(droplevels(ID))))
    }
    
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
    
    
    row <- data.frame(chr=chr,start= mean(tmp$start),end=mean(tmp$end),count=tmp$count[1], width=mean(tmp$width),n=length(levels(droplevels(tmp$ID))),perc=j,BLM=BLM,RECQL5=RECQL5,BLM_RECQL5=BLM_RECQL5,WT=WT)
    summary <- rbind(summary,row)
    
    if (printing){
      file.copy(files2transfer,readsdatapath)
      files2plot=paste0("DATA/rdata/",levels(tmp$ID),".RData")
      plotBreakpointsPerChr(files2plot,plotspath = datapath,chromosomes = c(chr))
    }
    if (cfs){
      file.copy(paste0("DATA/fragile_site_bed/",chr,"_fragile_site.bed"),datapath)
      
    }
  }
  return(summary)
}



################################################
        # Running Core Functions #
################################################

master <- function(printing = F,normalize="By_Library",export = F){
  breaks.all.files <- collectBreaksAllFiles(datapath="DATA/rdata/")
  
  hotspots <- breakpointHotspotter(breaks.all.files)
  
  summary <- savingAndPrinting(hotspots,hotpath="HOTSPOT_EVENTS/",printing = printing,normalize=normalize,export = export,genomeInstability = genomeInstability)

  message("\nI found ",nrow(summary), " hotspots.\n")
  message("\nThe average resolution is ", round(mean(summary$width),digits = -3),".\n")
  
  write.table(summary,"SUMMARY/hotspots.txt",col.names = T,row.names = F, quote = F,sep="\t")
  
  return(summary)
}


master(printing = T,normalize=F,export = T,genomeInstability=F)




  

