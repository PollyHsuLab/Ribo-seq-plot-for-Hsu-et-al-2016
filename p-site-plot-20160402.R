##################################################
# Only plot the riboseq reads in the CDS region  #
##################################################

#rm(list=ls())
setwd("~/Desktop/ORF_plots/")
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)

##################################################
#This file contains multiple functions to plot the gene models and the riboseq frame information.
#This gene model plot shows each exon with different colors to account for the frame shift in respect to the start position of the isofome.
#The riboseq counts also account for the frame shift in respect to the start position of the isofome (default is isoform 1 of the gene).

#Read in gtf (has to be gtf, no gff)
txdb <- makeTxDbFromGFF("~/Desktop/ForAraport/Araport11_20180126.gtf",
                        format="gtf",
                        dataSource="TAIR",
                        organism="Arabidopsis thaliana")

exonsByTx <- exonsBy(txdb,by='tx',use.names=T)
exonsByGene <- exonsBy(txdb,by='gene')
exonsByTxGene <- exonsBy(txdb,by=c('tx','gene'),use.names=T)
txByGene <- transcriptsBy(txdb,by='gene')
cdsByTx <-cdsBy(txdb, by="tx",use.names=T)
fiveUTR <- fiveUTRsByTranscript(txdb,use.names=T)
threeUTR <- threeUTRsByTranscript(txdb,use.names=T)
#Set a database with all cds information split by gene and transcripts
cds <- cdsBy(txdb, by=c("tx","gene"),use.names=TRUE)

############################
# Read in Riboseq results  #
############################
#The Root_P_sites_all.bed is derive from the p-site-all file of RiboTaper output
# you can just sum up the read count for each position of the p_sites_all file (in linux) to make the Root_P_sites_all.bed file

a <- read.delim("/Volumes/backup/Root_P_sites_all.bed",head=F,stringsAsFactors=F,sep="\t")#root data

head(a)
#   V1   V2   V3 V4 V5 V6
#   1  1 5694 5695  .  1  -
#   2  1 6918 6919  . 19  -
#   3  1 7044 7045  .  7  -
#   4  1 7156 7157  . 30  -
#   5  1 7159 7160  . 67  -
#   6  1 7162 7163  .  6  -

colnames(a) <- c("chr","position1","position","X","count","strand")
#remember tha the position1 and X are not required, you can discard that when you make the Root_P_sites_all.bed

ribo <- a[,c("count","chr","position","strand")]

############################
# Read in RNAseq results   #
############################
AtRNAseq <- "/Volumes/backup/RNA_R1+R2+R3.bam"
# This version of code load the bam file to the RAM so its slow and need to make sure that you have enough memory
# x is a class of GAlignments, Each row is a read
x <- readGAlignments(AtRNAseq)
#Coverage of the reads : this will generate a RleList
#Rle is a run-length encoding
xcov <- coverage(x)

############################
#plotRanges is a function to plot GRanges 
#p_site_plot_p plot periodicity of the p-site with respect to the transcripts
#p_site_plot_p2 is for plotting the uORFs, but you would need the uORF ranges in the gtf (which we did manually for the PNAS paper)
#plotGeneModel combines both plotRanges and p_site_plot_all functions
#firstInFramePSitePerExonPositive finds out the frame information with respect to isoform CDS start site when gene is on + strand
#firstInFramePSitePerExonNegative finds out the frame information with respect to isoform CDS start site when gene is on - strand

############################################
# plotRanges is a function to plot GRanges #
############################################
#x is the name of a gene model
plotRanges <- function(x,gene,uORF=NULL,shortest3UTR,ybottom = 0,main = deparse(substitute(x)),colCDS = "black",col3="white",col5="lightgrey",...) {
  
  if(x %in% names(cdsByTx)) {
    height <- 0.1
    xlim=ranges(unlist(exonsByTx[x]))
    xlimCds=ranges(unlist(cdsByTx[x]))
    
    if (x %in% names(fiveUTR)) {
      xlim5=ranges(unlist(fiveUTR[x]))
      rect(start(xlim5), ybottom, end(xlim5), ybottom + height, col = col5, border = "black",...)
    }
    if (length(unlist(exonsByTx[x]))>1) {
      GAPS <- gaps(unlist(exonsByTx[x]),start=NA)
      segments(x0 = start(GAPS),
               y0 = ybottom+height/2,
               x1 = start(GAPS)+width(ranges(GAPS))/2,
               y1 = ybottom+height,
               col = "black",lwd=1)
      segments(x0 = start(GAPS)+width(ranges(GAPS))/2,
               y0 = ybottom+height,
               x1 = end(GAPS),
               y1 = ybottom+height/2,
               col = "black",lwd=1)
    }
    #######################
    if (as.character(runValue(strand(exonsByTx[x])))=="+") {
      Frame <- firstInFramePSitePerExonPositive(x)

    }
    else {
      Frame <- firstInFramePSitePerExonNegative(x)

    }
    rect(start(xlimCds), ybottom, end(xlimCds), ybottom + height, col =c("black","black","black")[Frame] , border = "black",...)
    if (!is.null(uORF)) {
      uORF=paste0(uORF,".1","")
      xlim_uORF=ranges(unlist(cdsByTx[uORF]))
      rect(start(xlim_uORF), ybottom, end(xlim_uORF), ybottom + height, col ="yellow" , border = "black",...)
      }
    
    ########################
    
    if (x %in% names(threeUTR)) {
      xlim3=ranges(sort(unlist(threeUTR[x])))
      Length=length(unlist(threeUTR[x]))
      
      if (shortest3UTR <=50) {
        z=shortest3UTR
      }
      else{
        z=shortest3UTR/3
      }
      if (as.character(runValue(strand(exonsByTx[x])))=="+") {
        if (length(unlist(threeUTR[x]))==1) {
          polygon(x=c(start(xlim3), end(xlim3)-z,end(xlim3),end(xlim3)-z,start(xlim3)),y=c(ybottom+height,ybottom+height,ybottom+height/2,ybottom,ybottom), col = col3, border = "black",...)
        }
        else {
          rect(start(xlim3[1:Length-1]), ybottom, end(xlim3[1:Length-1]), ybottom + height, col = col3, border = "black",...)
          polygon(x=c(start(xlim3[Length]), end(xlim3[Length])-z,end(xlim3[Length]),end(xlim3[Length])-z,start(xlim3[Length])),y=c(ybottom+height,ybottom+height,ybottom+height/2,ybottom,ybottom), col = col3, border = "black",...)
        }
      }
      if (as.character(runValue(strand(exonsByTx[x])))=="-") {

        if (length(unlist(threeUTR[x]))==1) {
          polygon(x=c(start(xlim3),start(xlim3)+z,end(xlim3),end(xlim3),start(xlim3)+z),y=c(ybottom+height/2,ybottom+height,ybottom+height,ybottom,ybottom), col = col3, border = "black",...)
        }
        else {
          rect(start(xlim3[2:Length]), ybottom, end(xlim3[2:Length]), ybottom+height, col = col3, border = "black",...)
          polygon(x=c(start(xlim3[1]),start(xlim3[1])+z,end(xlim3[1]),end(xlim3[1]),start(xlim3[1])+z),y=c(ybottom+height/2,ybottom+height,ybottom+height,ybottom,ybottom), col = col3, border = "black",...)
        }
      }
    }
    axis(1)
  } 
  else {
    stop("Input transcript is not a coding gene in gtf/gff file.")
  }
}


###########################################################################
# p_site_plot_p function to plot periodicity according to transcript info #
###########################################################################


p_site_plot_p <- function(GeneName,isoform=1,CDSonly=FALSE) {
  #CDSonly=T, then only plot the reads in the CDS
  if(paste0(GeneName,".",isoform,sep = "") %in% names(cdsByTx)) {
    CDS <- cds[names(cds)==paste(GeneName,".",isoform,sep = ""),]
    #find ranges of exons 
    Exon <- exonsByGene[names(exonsByGene)==GeneName,]
    #Extract chromosome number from CDS object 
    chr=as.numeric(seqnames(unlist(CDS)))[1]
    #Extract strand information from CDS object 
    txStrand=as.character(strand(unlist(CDS)))[1]
    #Extract the CDS ranges
    cdsRanges = cdsByTx[names(cdsByTx)==paste0(GeneName,".",isoform,sep = ""),]
    #Extract most left position from the Exon object 
    txLeft <<-min(start(ranges(unlist(Exon))))
    #Extract most right position from the Exon object 
    txRight <<-max(end(ranges(unlist(Exon))))
    #Extract most left position from the CDS object 
    cdsLeft=min(start(ranges(unlist(CDS))))
    #Extract most right position from the CDS object 
    cdsRight=max(end(ranges(unlist(CDS))))
    ##Extract start site from CDS object 
    cdsStart=ifelse(txStrand=="+",as.numeric(min(start(ranges(CDS)))),as.numeric(max(end(ranges(CDS)))))
    cdsEnd=ifelse(txStrand=="+",as.numeric(max(end(ranges(CDS)))),as.numeric(min(start(ranges(CDS)))))
    
    s=data.frame()
    #Generate the sequences of positions in the transcript
    if(txStrand=="+") {
      sposition = sort(unlist(mapply(seq,start(unlist(cdsRanges)),end(unlist(cdsRanges)))),decreasing=F)
      sseq=seq(1,length(sposition),1)
      s1 = sposition[which(sseq%%3==1)]
      s2 = sposition[which(sseq%%3==2)]
      s3 = sposition[which(sseq%%3==0)]
    }
    else {
      sposition = sort(unlist(mapply(seq,start(unlist(cdsRanges)),end(unlist(cdsRanges)))),decreasing=T)
      sseq=seq(1,length(sposition),1)
      s1 = sposition[which(sseq%%3==1)]
      s2 = sposition[which(sseq%%3==2)]
      s3 = sposition[which(sseq%%3==0)]
    }
    #Extract riboseq reads in the region of the transcript 
    if (CDSonly==TRUE) {
      a <- ribo[ribo[,2]==chr & ribo[,3] > cdsLeft & ribo[,3] < cdsRight & ribo$strand==txStrand,]
    }
    else if(CDSonly==FALSE) {
      a <- ribo[ribo[,2]==chr & ribo[,3] > txLeft & ribo[,3] < txRight & ribo$strand==txStrand,]
    }
    
    a$frame <- as.factor(ifelse(a$position%in%s1,0,ifelse(a$position%in%s2,1,ifelse(a$position%in%s3,2,4))))
    YLIM <<- c(0,max(c(0,a$count)))

    plot(x=a$position,y=a$count,type="h",ylab="Count",
         xlim=c(txLeft,txRight),ylim=YLIM, #c(0,max(c(0,30))),
         col=c("red","#3366FF","#009900","darkgrey")[a$frame],
         lwd=1,xaxt = "n")
    axis(side=1, labels=FALSE, tck = -0.01)
    abline(v=cdsStart,lty=2,lwd=1)
    abline(v=cdsEnd,lty=2,lwd=1, col="darkgrey")
  }
  else {
    stop("Input transcript is not a coding gene in gtf/gff file.")
  }
}

p_site_plot_p2 <- function(uORF,isoform=1,CDSonly=FALSE) {
  #CDSonly=T, then only plot the reads in the CDS
  if(paste0(uORF,".",isoform,sep = "") %in% names(cdsByTx)) {
    CDS <- cds[names(cds)==paste(uORF,".",isoform,sep = ""),]
    #find ranges of exons 
    Exon <- exonsByGene[names(exonsByGene)==uORF,]
    #Extract chromosome number from CDS object 
    chr=as.numeric(seqnames(unlist(CDS)))[1]
    #Extract strand information from CDS object 
    txStrand=as.character(strand(unlist(CDS)))[1]
    #Extract the CDS ranges
    cdsRanges = cdsByTx[names(cdsByTx)==paste0(uORF,".",isoform,sep = ""),]
    #Extract most left position from the Exon object 
    #txLeft=min(start(ranges(unlist(Exon))))
    #Extract most right position from the Exon object 
    #txRight=max(end(ranges(unlist(Exon))))
    #Extract most left position from the CDS object 
    cdsLeft=min(start(ranges(unlist(CDS))))
    #Extract most right position from the CDS object 
    cdsRight=max(end(ranges(unlist(CDS))))
    ##Extract start site from CDS object 
    cdsStart=ifelse(txStrand=="+",as.numeric(min(start(ranges(CDS)))),as.numeric(max(end(ranges(CDS)))))
    cdsEnd=ifelse(txStrand=="+",as.numeric(max(end(ranges(CDS)))),as.numeric(min(start(ranges(CDS)))))
    
    s=data.frame()
    #Generate the sequences of positions in the transcript
    if(txStrand=="+") {
      sposition = sort(unlist(mapply(seq,start(unlist(cdsRanges)),end(unlist(cdsRanges)))),decreasing=F)
      sseq=seq(1,length(sposition),1)
      s1 = sposition[which(sseq%%3==1)]
      s2 = sposition[which(sseq%%3==2)]
      s3 = sposition[which(sseq%%3==0)]
    }
    else {
      sposition = sort(unlist(mapply(seq,start(unlist(cdsRanges)),end(unlist(cdsRanges)))),decreasing=T)
      sseq=seq(1,length(sposition),1)
      s1 = sposition[which(sseq%%3==1)]
      s2 = sposition[which(sseq%%3==2)]
      s3 = sposition[which(sseq%%3==0)]
    }
    #Extract riboseq reads in the region of the transcript 
    if (CDSonly==TRUE) {
      a <- ribo[ribo[,2]==chr & ribo[,3] > cdsLeft & ribo[,3] < cdsRight & ribo$strand==txStrand,]
    }
    else if(CDSonly==FALSE) {
      a <- ribo[ribo[,2]==chr & ribo[,3] > txLeft & ribo[,3] < txRight & ribo$strand==txStrand,]
    }
    
    a$frame <- as.factor(ifelse(a$position%in%s1,0,ifelse(a$position%in%s2,1,ifelse(a$position%in%s3,2,4))))
    
    plot(xlim=c(txLeft,txRight),YLIM, xaxt = "n",yaxt = "n",type="n")
    polygon(x=c(cdsStart,cdsEnd,cdsEnd,cdsStart),y=c(0,0,max(YLIM),max(YLIM)),border = "white",col="white")
    lines(x=a$position,y=a$count,type="h",ylab="Count",xlim=c(txLeft,txRight),YLIM,col=c("red","#3366FF","#009900","dark grey")[a$frame],lwd=1,xaxt = "n")
    abline(v=cdsStart,lty=2,lwd=1)
    abline(v=cdsEnd,lty=2,lwd=1, col="darkgrey")
  }
  else {
    stop("Input transcript is not a coding gene in gtf/gff file.")
  }
}

plotGeneModel <- function(gene,uORF){
  isoforms <- length(unlist(txByGene[gene]))
  generanges <- ranges(unlist(exonsByGene[gene]))
  genelim <- c(min(start(generanges)), max(end(generanges)))
  GeneIsoformNumber <- length(unlist(txByGene[gene]))

  plot.new()
  yAxis <- (isoforms*0.2+0.1)
  plot.window(genelim,c(0,yAxis))
  for (i in 1:isoforms) {
    x=paste0(gene,".",i,"")
    if (x %in% names(threeUTR)) {
      if (i==1) {
        shortest3UTR <- min(sapply(1:isoforms, function(i) sum(width(unlist(threeUTR[paste0(gene,".",i,"")])))))
        plotRanges(x=paste0(gene,".",i,""),gene,uORF,shortest3UTR,ybottom=(yAxis-0.2*i)) 
        
      }
      else {
        shortest3UTR <- min(sapply(1:isoforms, function(i) sum(width(unlist(threeUTR[paste0(gene,".",i,"")])))))
        plotRanges(x=paste0(gene,".",i,""),gene,uORF,shortest3UTR,yaxt="n",ybottom=(yAxis-0.2*i)) 
        
      }}
                  
    else {
      plotRanges(x=paste0(gene,".",i,""),gene,uORF,ybottom=(yAxis-0.2*i))
    }
  }
  
}

#######################################################
firstInFramePSitePerExonPositive <- function(x){
  nCDS =length(unlist(cdsByTx[x]))
  YFGrange = as.numeric(unlist(ranges(unlist(cdsByTx[x]))))
  Range = seq(1,length(YFGrange))
  Seq3 = as.factor((Range-1)%%3)
  dfs = data.frame(YFGrange,Range,Seq3)
  dfs3 = data.frame()
  Listdf = split(dfs,rep(1:nCDS,width(unlist(cdsByTx[x]))))
  for (i in 1:nCDS) {
    dfs3 = rbind(dfs3,Listdf[[i]][Listdf[[i]]$Seq3==0,][1,])
  }
  dfs3$frame <- as.factor((dfs3$YFGrange- dfs3$YFGrange[1])%%3)
  levels(dfs3$frame)=c(0,1,2)
  return(dfs3$frame)
}

firstInFramePSitePerExonNegative <- function(x){
  nCDS =length(unlist(cdsByTx[x]))
  YFGrange = sort(as.numeric(unlist(ranges(unlist(cdsByTx[x])))),decreasing=T)
  Range = seq(1,length(YFGrange))
  Seq3 = as.factor((Range-1)%%3)
  dfs = data.frame(YFGrange,Range,Seq3)
  Listdf = split(dfs,rep(1:nCDS,width(unlist(cdsByTx[x]))))
  dfs3 = data.frame()
  for (i in 1:nCDS) {
    dfs3 = rbind(dfs3,Listdf[[i]][Listdf[[i]]$Seq3==0,][1,])
  }
  dfs3$frame <- as.factor((dfs3$YFGrange[1]-dfs3$YFGrange)%%3)
  levels(dfs3$frame)=c(0,1,2)
  return(dfs3$frame)
}

###############################################

PLOT <-function(YFG,uORF=NULL,CDSonly=FALSE,isoform=1) {
  par(mfrow=c(3,1),mar=c(0.2,0.2,0.2,0.2),oma=c(3,2,4,2))
  chr <- as.numeric(substr(YFG,3,3))
  generanges <- ranges(unlist(exonsByGene[YFG]))
  GR <<- GRanges(seqnames=chr,IRanges(min(start(generanges)), max(end(generanges))),strand=strand(unlist(exonsByGene[YFG]))[1])
  Gtx <<- as.numeric(xcov[[chr]][ranges(GR)])
  isoforms <- length(unlist(txByGene[YFG]))
  layout(matrix(c(1,1,2,2,3,3),3,2,byrow=TRUE), widths=c(6,6,6), heights=c(1.5,1.5,0.35*isoforms))
  plot(Gtx,type="h",col="grey",lwd=1,xaxt='n',ylim=c(0,max(Gtx)+2))
  lines(x=c(1,length(Gtx)),y=c(0,0),col="white",lwd=2)
  legend("topleft",legend="RNAseq",bty="n",cex=1.5)
  p_site_plot_p(YFG,CDSonly=CDSonly,isoform=isoform)
  par(new=TRUE)
  if (!is.null(uORF)) {p_site_plot_p2(uORF,CDSonly=TRUE,isoform=1)}
  legend("topleft",legend="Riboseq",bty="n",cex=1.5)
  plotGeneModel(YFG,uORF)
  mtext(YFG,NORTH<-3,line=0.4, cex=1.2, col="black", outer=TRUE)
}


PLOT(YFG ="AT2G46980",isoform=3)
PLOT(YFG ="AT2G47000",isoform=1) #Good example
PLOT(YFG ="AT2G47060",isoform=1)
PLOT(YFG ="AT2G47070",isoform=1)

