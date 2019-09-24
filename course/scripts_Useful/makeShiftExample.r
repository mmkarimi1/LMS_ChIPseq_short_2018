library(GenomicAlignments)
peaks <- ChIPQC:::GetGRanges("~/Downloads/bigchip/Peaks/MACS2/SRR568130/SRR568130_withInput_SRR502663_peaks.narrowPeak")
peaks <- peaks[order(peaks$X61.48394,decreasing=T),]
readsP1 <- as(readGAlignments("~/Downloads/bigchip/BAMs/Sorted_SRR568130.bam",
                         param=ScanBamParam(which=peaks[1,])),"GRanges")

minStart <- min(start(readsP1))-1
start(readsP1) <- start(readsP1)-minStart
end(readsP1) <- end(readsP1)-minStart
readsP1_GA <- as(readsP1,"GAlignments")

library(Gviz)
mycch12ReadsAT_Pos <- AlignmentsTrack(resize(readsP1[strand(readsP1) == "+"],1))
mycch12ReadsAT_Neg <- AlignmentsTrack(resize(readsP1[strand(readsP1) == "-"],1))

library(gridExtra)
myCor <- vector()
for(x in 1:300){

mycch12ReadsAT_Pos <- AlignmentsTrack(shift(readsP1[strand(readsP1) == "+"],x),fill.coverage="Darkred")
mycch12ReadsAT_Neg <- AlignmentsTrack(readsP1[strand(readsP1) == "-"],,fill.coverage="Darkgreen")
myCor[x] <- cor(coverage(shift(readsP1[strand(readsP1) == "+"],x))[[2]],coverage(readsP1[strand(readsP1) == "-"])[[2]])
png(paste0("test",x,".png"),width=2000,height=1000)
plotTracks(c(mycch12ReadsAT_Pos,mycch12ReadsAT_Neg),type="coverage",from=0,to=1000)
dev.off()
png(paste0("testCor_new",x,".png"),width=2000,height=1000)
p <- ggplot(data.frame(Shift=1:x,Correlation=myCor),aes(x=Shift,y=Correlation))+geom_point(size=10,alpha=0.8)+xlim(c(0,400))+ylim(c(0,1))+stat_smooth(size=8,method="loess")+theme(axis.title=element_text(size = 40),axis.title.y=element_text(angle=0))
print(p)
dev.off()
}
