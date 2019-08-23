
library("tidyverse")
library("reshape2")
library("DESeq2")
library("gplots")
library("fgsea")

# read counts
tmp<-read.table("3col.tsv",header=F)
x<-as.matrix(acast(tmp, V2~V1, value.var="V3", fun.aggregate = sum))
x<-as.data.frame(x)
accession<-sapply((strsplit(rownames(x),"\\|")),"[[",2)
symbol<-sapply((strsplit(rownames(x),"\\|")),"[[",6)
x$geneid<-paste(accession,symbol)

xx<-aggregate(. ~ geneid,x,sum)
rownames(xx)<-xx$geneid
xx$geneid=NULL
xx<-round(xx)

pdf("hmec1_MDSplot.pdf")
plot(cmdscale(dist(t(xx))), xlab="Coordinate 1", ylab="Coordinate 2", type = "n") 
text(cmdscale(dist(t(xx))), labels=colnames(xx)) 
dev.off()
# curate the samplesheet
samplesheet<-as.data.frame(colnames(xx))
samplesheet$hcy<-c(0,0,0,1,1,1)
rownames(samplesheet)<-samplesheet[,1]
samplesheet[,1]=NULL

####################################################
# ctrl vs hcy
####################################################
y<-xx[which(rowSums(xx)/ncol(xx)>=(20)),]

dds <- DESeqDataSetFromMatrix(countData = y , colData = samplesheet, design = ~ hcy )
res <- DESeq(dds)
z<- results(res)
vsd <- vst(dds, blind=FALSE)
zz<-cbind(as.data.frame(z),assay(vsd))
dge<-as.data.frame(zz[order(zz$pvalue),])
write.table(dge,file="hmec1_hcy_deseq.tsv",quote=F,sep="\t")
rnk<-as.data.frame( sign(dge$log2FoldChange) * (-log(dge$pvalue + 1E-307) ))
rownames(rnk)<-row.names(dge)
colnames(rnk)="Score"
write.table(rnk,file="hmec1_hcy_deseq.rnk",sep='\t',quote=F)

#some plots
pdf("hmec1_hcy_plots.pdf")
sig<-subset(dge,padj<0.05)
SIG=nrow(sig)
DN=nrow(subset(sig,log2FoldChange<0))
UP=nrow(subset(sig,log2FoldChange>0))
HEADER=paste("Ctrl vs 200 uM Hcy:", SIG , "DGEs,", UP ,"upregulated,", DN, "downregulated")

plot(log2(dge$baseMean),dge$log2FoldChange,cex=0.5, xlab="log2 base mean", 
 ylim=c(-3,3),ylab="log2 fold change" ,pch=19,col="#838383")

points(log2(sig$baseMean),sig$log2FoldChange,cex=0.5,pch=19,col="red")
mtext(HEADER)
top<-head(sig,20)
#text(log2(top$baseMean)+1, top$log2FoldChange, labels = rownames(top),cex=0.7)
#volcano plot
plot(dge$log2FoldChange, -log2(dge$pvalue) ,cex=0.5, xlim=c(-3,3),xlab="log2 fold change", ylab="-log2 p-value" ,pch=19,col="#838383")
points(sig$log2FoldChange, -log2(sig$pvalue),cex=0.5,pch=19,col="red")
#text(top$log2FoldChange+0.5, -log2(top$pvalue), labels = rownames(top),cex=0.7)
mtext(HEADER)
# top N gene heatmap
colfunc <- colorRampPalette(c("blue", "white", "red"))
heatmap.2(  as.matrix(dge[1:50,c(7:ncol(dge))]), col=colfunc(25),scale="row",
 trace="none",margins = c(6,20), cexRow=.6, cexCol=.6,  main="Top 50 genes")
dev.off()

####################################################
# FGSEA
####################################################

dge$SYMBOL<-sapply((strsplit(rownames(dge)," ")),"[[",2)

res2 <- dge %>% 
  dplyr::select(SYMBOL,stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))
res2

stat<- deframe(res2)

download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip",
 destfile="ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")

gsets<-gmtPathways("ReactomePathways.gmt")

fgseaRes <- fgsea(pathways=gsets, stats=stat, minSize=10, nperm=100000)

fgseaRes<-fgseaRes[order(fgseaRes$pval),]


# peek at upregulated pathways
head(subset(fgseaRes,ES>0),10)

# peek at downregulated pathways
head(subset(fgseaRes,ES<0),10)

pdf("hmec1_pathway_charts.pdf")

psig<-subset(fgseaRes,padj<=0.05)
plot(fgseaRes$ES,-log10(fgseaRes$pval),pch=19,cex=0.8,xlab="ES",ylab="-log10(p-value)")
points(psig$ES,-log10(psig$pval),pch=19,cex=0.8,xlab="ES",ylab="-log10(p-value)",col="red")
TOTAL=nrow(fgseaRes)
SIG=nrow(psig)
UP=length(which(psig$ES>0))
DN=length(which(psig$ES<0))
HEADER=paste(TOTAL,"gene sets examined,",SIG,"FDR<0.05,",UP,"up-regulated,",DN,"down-regulated")
mtext(HEADER)

# barchart top effect size FDR<0.05
par(mfrow = c(2, 1)) 
sig_up<-subset(psig,ES>0)
sig_dn<-subset(psig,ES<0)
sig_up<-head(sig_up[order(-sig_up$ES),],20)
sig_dn<-head(sig_dn[order(-sig_dn$ES),],20)

par( mar=c(5,28,4,2)) ; barplot(rev(sig_up$ES),horiz=T,names.arg=rev(sig_up$pathway),las=2, cex.names=0.6 , cex.axis=0.6, xlab="ES",main="up-regulated sets",cex.main=0.6)
par( mar=c(5,28,4,2)) ; barplot(rev(sig_dn$ES),horiz=T,names.arg=rev(sig_dn$pathway),las=2, cex.names=0.6, cex.axis=0.6, xlab="ES",main="down-regulated sets",cex.main=0.6)
dev.off()			`
# barchart top significant


# Pathway Table
psig$leadingEdge <- vapply(psig$leadingEdge, paste, collapse = ", ", character(1L))
write.table(psig, file = "Hmec1_pathway_tables.csv", sep = ",")
