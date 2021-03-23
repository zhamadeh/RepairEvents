################################################
            # Packages #
################################################
library(tidyverse)
library(GenomicRanges)
library(scales)

################################################
          # Plotting #
################################################

summary = read.table("SUMMARY/hotspots.txt",header=T)
densitySummary = read.table("SUMMARY/densityPvalueSummary.txt",header=T)

################################################
#PLOTTING WITH TRANSPARENT BACKGROUND
################################################

transparentBackground <- function(p,filename){
  p <- p +
    theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    ) 
  ggsave(p, filename = paste0("PLOTS/",filename),  bg = "transparent")
}


################################################
# DENSITY SUMMARY FOR ALL CHROMOSOMES AND CHR 1
################################################

max_xy <-max(densitySummary$ds)
densitySummary$mp_Mb = densitySummary$mp/1e+06
densitySummary$chr <- factor(densitySummary$chr,levels = c("chr1" , "chr2",  "chr3",  "chr4" , "chr5" , "chr6",  "chr7",  "chr8" , "chr9", "chr10" ,"chr11" ,"chr12" ,"chr13", "chr14", "chr15" ,"chr16" ,"chr17", "chr18","chr19" ,"chr20", "chr21" ,"chr22", "chrX" ))

densityALL = ggplot(densitySummary)+geom_point(aes(mp_Mb,ds,color=type),size=1)+facet_wrap(~chr,scales="free")+
  scale_y_continuous(limits = c(0,max_xy)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.title =element_text(size=28),
        text = element_text(size=20),
        legend.position = c(0.7, 0.07))+
  labs(y="DENSITY",x="CHROMOSOME POSITION")+
  guides(color = guide_legend(override.aes = list(size=10)))

densityChr1 = ggplot(dplyr::filter(densitySummary,chr=="chr1"))+geom_line(aes(mp_Mb,ds,color=type),size=4)+facet_wrap(~chr,scales="free")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=18),
        axis.title =element_text(size=28),
        text = element_text(size=20),
        legend.position = c(0.07, 0.7))+
  labs(y="DENSITY",x="CHROMOSOME POSITION")+
  guides(color = guide_legend(override.aes = list(size=8)))



################################################
# PLOTTING HOTSPOTS PER CHROMOSOME AND EVENTS PER HOTSPOT
################################################

summary$count<- seq(1:nrow(summary))
summarize <- gather(summary, gene,freq,c(BLM,RECQL5,BLM_RECQL5,WT))
summarize$chr<-factor(summarize$chr,levels=c("chr1"  ,"chr2" ,"chr3" , "chr4"  ,"chr5" , "chr7",  "chr8" , "chr9"  ,"chr10", "chr11" ,"chr12" ,"chr13", "chr14"
                                             ,"chr15" ,"chr16", "chr17", "chr18", "chr19" , "chr20", "chr21", "chr22", "chrX"))

summarize = summarize[order(summarize$chr),]
summarize$chr<-as.character(summarize$chr)
summarize$count<-as.character(summarize$count)

eventsPerHotspot = ggplot(summarize)+geom_col(aes(count,freq,fill=gene))+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  labs(y="Number of events/library",x="Hotspots")

hotspotsPerChromosome =ggplot(summarize)+geom_col(aes(chr,freq,fill=gene))+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  labs(y="# OF HOTSPOTS",x="CHROMOSOME")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=28),
        text = element_text(size=20),
        legend.position = "none")

summarize$chr <- as.factor(summarize$chr)





################################################
# COMMON FRAGILE SITE PLOTTING
################################################



cfs <-read.table("DATA/fragile_site_bed/CFS_human.bed",fill=T) %>% select(c(V1,V2,V3))
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

cfsHotspots <- ggplot()+ geom_bar(s,mapping=aes(chr,n,fill=type),stat="identity",position=  position_stack(reverse=T))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  labs(x="Chromosome",y="Count")







################################################
# INVERSION PLOTTING
################################################



summarizeToPlot = filter(summarize,freq>0)
summarizeToPlot$gene=as.factor(summarizeToPlot$gene)
numOfLibsPerGene$gene<-as.character(numOfLibsPerGene$gene)
numOfLibsPerGene[2,1]="BLM_RECQL5"
merge =merge(summarizeToPlot,numOfLibsPerGene,by="gene")
merge$n..=as.numeric(merge$n..)
merge$norm = merge$freq/merge$n..
merge$count <- as.factor(merge$count)

merge$count=as.numeric(merge$count)
merge = merge[order(merge$chr),]
merge$count <- as.factor(merge$count)
merge = filter(merge,count!="13" & count!= "27")
inversions = ggplot(merge)+geom_point(aes(gene,count,size=norm,color=norm))+
  labs(size="ALLELE FREQUENCY")+
  theme(legend.position="top",
        legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        axis.text.x.bottom  = element_text(size=20),
        axis.title.x=element_blank(),
        axis.title =element_text(size=25),
        text = element_text(size=18),
        axis.text.x=element_text(angle=60,hjust=0.5,vjust=0.5))+
  scale_size_continuous(range=c(1,7))+
  guides(color=F)+
  labs(y="INVERSIONS")


summarizeToPlot$gene=as.factor(summarizeToPlot$gene)

inversionWidth = ggplot(summarizeToPlot)+geom_point(aes(gene,width,size=freq,color=freq))+scale_y_continuous(trans = 'log10',
                                                                                            breaks = trans_breaks('log10', function(x) 10^x),
                                                                                            labels = trans_format('log10', math_format(10^.x)))+
  theme(legend.title = element_blank(),
        legend.position="top",
        axis.title.x=element_blank(),
        axis.title =element_text(size=18),
        text = element_text(size=18),
        axis.text.x=element_text(angle=60,hjust=0.5,vjust=0.5))+
  guides(color=F)+
  labs(y="INVERSION WIDTH")



