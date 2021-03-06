

savingAndPrinting(hotspots)



savingAndPrinting <- function(hotspots){
  
  hotpath="SCE_HOTSPOTS"
  if (!file.exists(hotpath) ) { dir.create(hotpath)}
  
  hotspots$count=as.factor(hotspots$count)
  hotspots$ID <- as.factor(hotspots$ID)
  
  summary <- data.frame(chr=c(),start= c(),end=c(), width=c(),n=c(),perc=c(),BLM=c(),RECQL5=c(),BLM_RECQL5=c(),WT=c())
  
  
  
  for (level in levels(hotspots$count)){
    tmp <- filter(hotspots, count==level)
    chr = as.character(tmp[1,1])
    
    tmp$ID <- droplevels(tmp$ID)
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
    message("Found a hotspot on ", chr," in ", j,"% of 317 cells.")
    
    
    if (perc[which.max(perc$perc),]$perc > 90){
      message("This inversion is unique to ", perc[which.max(perc$perc),]$gene, " making up ",perc[which.max(perc$perc),]$perc, " of the libraries")
    }
    
    
    datapath=paste0(hotpath,"/",chr,"-",j,"/")
    if (!file.exists(datapath) ) { dir.create(datapath)}
    readsdatapath=paste0(hotpath,"/",chr,"-",j,"/reads/")
    if (!file.exists(readsdatapath) ) { dir.create(readsdatapath)}
    
    bed=tmp
    export(bed,paste0(datapath,"breakpoints.bed"),format = "gff3")
    write.table(as.data.frame(perc),paste0(datapath,"genotype.txt"),row.names = F,col.names = T,quote = F,sep="\t")
    #export(bed,paste0(datapath,"breakpoints.bed"),format = "bed")
    
    if (length(tmp$ID)>25){
      files2transfer=paste0("DATA/browserfiles/",tmp$ID[1:25],"_reads.bed.gz")
    } else {
      files2transfer=paste0("DATA/browserfiles/",tmp$ID,"_reads.bed.gz")
    }
    file.copy(files2transfer,readsdatapath)
    
    if ("BLM" %in% perc$gene){
      BLM=as.numeric(perc[perc$gene=="BLM",2])
    } else { BLM=0 }
    if ("RECQL5" %in% perc$gene){
      RECQL5=as.numeric(perc[perc$gene=="RECQL5",2])
    } else { RECQL5=0 }
    if ("BLM/RECQL5" %in% perc$gene){
      BLM_RECQL5=as.numeric(perc[perc$gene=="BLM/RECQL5",2])
    } else { BLM_RECQL5=0 }
    if ("WT" %in% perc$gene){
      WT=as.numeric(perc[perc$gene=="WT",2])
    } else { WT=0 }
    
    row <- data.frame(chr=chr,start= mean(tmp$start),end=mean(tmp$end), width=mean(tmp$width),n=nrow(tmp),perc=j,BLM=BLM,RECQL5=RECQL5,BLM_RECQL5=BLM_RECQL5,WT=WT)
    summary <- rbind(summary,row)
    
    
    files2plot=paste0("DATA/rdata/",levels(tmp$ID),".RData")
    plotBreakpointsPerChr(files2plot,plotspath = datapath,chromosomes = c(chr))
  }
}
