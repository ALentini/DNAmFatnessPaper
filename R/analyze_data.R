options(stringsAsFactors = F)

###################
## Load raw data ##
###################

library(minfi)
targets <- read.metharray.sheet("Data")
RGSet <- read.metharray.exp(targets=targets)
# Exclude outlier sample
RGSet2 <- RGSet[,!grepl("100896160012_R04C01",targets$Basename)]

###################
## Load metadata ##
###################

library(data.table)
meta <- fread("metadata.tsv")
# Order rows by 450k data
meta <- meta[match(pData(RGSet2)$Sample_Name, ID)]
# Probe annotation from https://doi.org/10.1093/nar/gkw967
hm450k <- fread("hm450.hg19.manifest.tsv.gz")

######################
## Preproccess data ##
######################

## Normalize data using functional normalization
GRset.funnorm <- preprocessFunnorm(RGSet2)
# Filter out bad probes
GRset.filter <- GRset.funnorm[!hm450k[match(featureNames(GRset.funnorm),probeID)]$MASK_general,]
# Restrict to autosomes
GRset.auto <- GRset.filter[hm450k[match(featureNames(GRset.filter),probeID)]$CpG_chrm %in% paste0("chr",1:22),]
GRset.anno <- as.data.table(getAnnotation(GRset.auto))

get.probe <- function(x){
  idx <- sapply(strsplit(GRset.anno$UCSC_RefGene_Name,";"),function(y) any(y == x) )
  GRset.anno$Name[idx]
}

b <- getBeta(GRset.auto)
m <- getM(GRset.auto)

##########################
## Plot global profiles ##
##########################

## HClustering
library(ggplot2)
library(ggdendro)
library(cowplot)

b_cl <- getBeta(GRset.auto)
colnames(b_cl) <- meta$ID

hc <- hclust(dist(t(b_cl)))

#Theme for plotting continuous heatmaps
p_cont <- ggplot(meta, aes(x=factor(ID),y=1)) + scale_x_discrete(limits=colnames(b_cl)[hc$order]) + scale_fill_gradient(low="white",high="red") + theme_void() + theme(legend.position = "top", axis.title.x = element_blank(),axis.title.y = element_blank())

p_cl <- list(
  #Dendrogram
  p_cl = ggdendrogram(hc),
  #Slide
  p1 = ggplot(as.data.frame(pData(GRset.auto)), aes(x=factor(Sample_Name),y=1,fill=factor(Slide))) + scale_x_discrete(limits=colnames(b_cl)[hc$order]) + geom_tile() + theme_void() + theme(legend.position = "top", axis.title.x = element_blank(),axis.title.y = element_blank()),
  #Sex
  p2 = ggplot(meta, aes(x=factor(ID),y=1,fill=factor(Sex_0_F))) + scale_x_discrete(limits=colnames(b_cl)[hc$order]) + geom_tile() + theme_void() + theme(legend.position = "top", axis.title.x = element_blank(),axis.title.y = element_blank()),
  #BMI_1v
  p3 = p_cont + geom_tile(aes(fill=BMI_1w)),
  #FMI_1v
  p4 = p_cont + geom_tile(aes(fill=FMI_1w)),
  #FFMI_1v
  p5 = p_cont + geom_tile(aes(fill=FFMI_1w)),
  #BMI_PA
  p6 <- p_cont + geom_tile(aes(fill=BMI_Pat)),
  #FMI_PA
  p7 = p_cont + geom_tile(aes(fill=FMI_Pat)),
  #FFMI_PA
  p8 = p_cont + geom_tile(aes(fill=FFMI_Pat)),
  #BMI_v32_MA
  p9 = p_cont + geom_tile(aes(fill=BMI_w32_Mat)),
  #FMI_MA
  p10 = p_cont + geom_tile(aes(fill=FMI_Mat)),
  #FFMI_MA
  p11 = p_cont + geom_tile(aes(fill=FFMI_Mat))
)

ggsave("plots/HClust.pdf",width = 8,height = 12,
  plot_grid(plotlist = p_cl, align = "v", ncol = 1, rel_heights = c(1/2, rep(1/20,2),rep(1/14,9)))
)

## PCA plot
library(ggfortify)
pca <- prcomp(t(b))

ggsave("plots/PCA_BMI.pdf",width = 12,height = 12,
  gridExtra::grid.arrange(
    autoplot(pca,data=meta,size="BMI_1w",colour="Sex_0_F") + theme_bw(),
    autoplot(pca,data=meta,size="FMI_1w",colour="Sex_0_F") + theme_bw(),
    autoplot(pca,data=meta,size="FFMI_1w",colour="Sex_0_F") + theme_bw(),
    autoplot(pca,data=meta,colour="Sex_0_F",size='BMI_prepreg') + theme_bw(),
    autoplot(pca,data=meta,colour="Sex_0_F",size='BMI_w32_Mat') + theme_bw(),
    autoplot(pca,data=meta,colour="Sex_0_F",size='BMI_Pat') + theme_bw()
  )
)

###################
## Identify DMPs ##
###################

get.dmp <- function(x,pheno=character(), type, ...){
  require(minfi)
  dmp <- dmpFinder(x,pheno=meta[[pheno]], type = type)
  return(dmp)
}

# Type of each variable
type <- setNames(c("categorical", rep("continuous",9) ),c("Low_vs_high_infant_body_fatness","BMI_1w","FMI_1w","FFMI_1w","BMI_Pat","FMI_Pat","FFMI_Pat","BMI_w32_Mat","FMI_Mat","FFMI_Mat"))

dmp <- lapply(c("Low_vs_high_infant_body_fatness","BMI_1w","FMI_1w","FFMI_1w","BMI_Pat","FMI_Pat","FFMI_Pat","BMI_w32_Mat","FMI_Mat","FFMI_Mat") ,function(x) get.dmp(b,x,type[x]) )
names(dmp) <- c("Low_vs_high_infant_body_fatness","BMI_1w","FMI_1w","FFMI_1w","BMI_Pat","FMI_Pat","FFMI_Pat","BMI_w32_Mat","FMI_Mat","FFMI_Mat")

# Not normally distributed
dmp_homa <- dmpFinder(b,log(meta$HOMA_IR),type="continuous")

#################
## Filter DMPs ##
#################

# Remove invariable probes to increase power
# Method from: https://doi.org/10.1186/s13148-017-0320-z
# However, more stringent cutoff of 20% variability ~95th percentile
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
b_var <- apply(b,1,Variation)
b_var2 <- row.names(b)[!b_var<0.2]

# Plot variability
pdf("plots/DMP_variability_cutoff.pdf",8,8)
plot(quantile(b_var,seq(0,1,0.01)),type="l",ylim=c(0,1),las=1, ylab="CpG variability (beta)", xlab= "Percentile");abline(h=0.2,col="red",lwd=2);text(x=10,y=0.2,label="Variability cutoff",pos = 3,col="red")
dev.off()

# Filter invariant probes and recalculate FDR
filter.dmp <- function(x){
  x_f <- x[row.names(x) %in% b_var2,]
  x_f$qval <- p.adjust(x_f$pval,"BH")
  return(x_f)
}
dmp_filter <- lapply(dmp,filter.dmp)
dmp_homa_filter <- filter.dmp(dmp_homa)

# Plot number of probes P<0.05
library(cowplot)
ggsave("plots/DMP_probes_signif.pdf",width = 4,height = 4,qplot(y=c(sapply(dmp_filter,function(x) sum(x$pval < 0.05) ),sum(dmp_homa_filter$pval < 0.05) ),x=c(names(dmp_filter),"HOMA_IR"), fill=c(gsub(".*_","",names(dmp_filter)),"HOMA_IR"),geom="col") + coord_flip() + theme_cowplot() + theme(axis.title.y = element_blank(), legend.position = "top", legend.title = element_blank()) + ylab("Probes P-value < 0.05"))

#dmp_melt <- melt(sapply(dmp_filter,function(x) x$pval ))
#ggplot(dmp_melt,aes(y=-log10(value),x=Var2)) + geom_violin() + geom_boxplot(width=0.1,outlier.stroke = NA) + coord_flip()

########################
## DMP manhattan plot ##
########################

library(ggman)
library(cowplot)
# Add positional information for manhattan plot
format.dmp <- function(x){
  df <- data.frame(cpg=row.names(x),pval=x$pval,qval=x$qval)
  df <- cbind(df,hm450k[match(row.names(x),probeID),c("CpG_chrm", "CpG_beg")])
  df$CpG_chrm <- gsub("chr","",df$CpG_chrm)
  return(df)
}
dmp_man <- lapply(dmp_filter, format.dmp )
dmp_homa_man <- format.dmp(dmp_homa_filter)

# Plot infant body fatness manhattan plot
ggsave("plots/Manhattan_fatness.pdf",width = 16,height = 4,
  gridExtra::grid.arrange(nrow=1,
    ggman(dmp_man[[1]],snp='cpg',bp='CpG_beg',chrom='CpG_chrm',pvalue='pval',pointSize = 2, sigLine = NA) + scale_colour_grey(start=0.5,end=0.6) + labs(title=names(dmp_man)[1],y="-log10 P-value") + coord_cartesian(ylim=c(0,4)) + geom_hline(yintercept = -log10(0.05),col="red",lwd=1) + theme_cowplot(10),
    ggman(dmp_man[[1]],snp='cpg',bp='CpG_beg',chrom='CpG_chrm',pvalue='qval',pointSize = 2, sigLine = NA) + scale_colour_grey(start=0.5,end=0.6) + labs(title=names(dmp_man)[1],y="-log10 FDR") + coord_cartesian(ylim=c(0,4)) + geom_hline(yintercept = -log10(0.05),col="red",lwd=1) + theme_cowplot(10)
  )
)

# Plot HOMA_IR manhattan plot
ggsave("plots/Manhattan_HOMA_IR.pdf",width = 16,height = 4,
gridExtra::grid.arrange(nrow=1,
  ggman(dmp_homa_man,snp='cpg',bp='CpG_beg',chrom='CpG_chrm',pvalue='pval',pointSize = 2, sigLine = NA) + scale_colour_grey(start=0.5,end=0.6) + labs(title="HOMA_IR",y="-log10 P-value") + coord_cartesian(ylim=c(0,4)) + geom_hline(yintercept = -log10(0.05),col="red",lwd=1) + theme_cowplot(10),
  ggman(dmp_homa_man,snp='cpg',bp='CpG_beg',chrom='CpG_chrm',pvalue='qval',pointSize = 2, sigLine = NA) + scale_colour_grey(start=0.5,end=0.6) + labs(title="HOMA_IR",y="-log10 FDR") + coord_cartesian(ylim=c(0,4)) + geom_hline(yintercept = -log10(0.05),col="red",lwd=1) + theme_cowplot(10)
))

# Plot parent variables
p_list <- lapply(names(dmp_man),function(x){
  require(ggman)
  df <- dmp_man[[x]]
  p <- ggman(df,snp='cpg',bp='CpG_beg',chrom='CpG_chrm',pvalue='pval',sigLine=-log10(0.05), pointSize=2)
  # If any value passes FDR<0.05, color the point red
  if(any(df$qval<0.05)){
    ls <- df$cpg[df$qval<0.05]
    p <- ggmanHighlight(p,highlight=ls,colour = "red")
  }else{
    p <- p + scale_colour_grey(start = 0.5,end = 0.6)
  }
  # Add labels to top 5 CpGs (not used)
  # p <- ggmanLabel(p, labelDfm = head(df,5),snp='cpg',label='cpg',type = "text", colour="black")
  p <- p + coord_cartesian(ylim=c(0,6)) + ylab("-log10 P-value") + ggtitle(x) + theme_cowplot(10)
  return(p)
})
# Plot in grid
library(cowplot)
ggsave("plots/Manhattan_pval.pdf",width = 24,height = 12,
  plot_grid(plotlist = p_list[-1])
)

# Plot FDR
q_list <- lapply(names(dmp_man),function(x){
  require(ggman)
  df <- dmp_man[[x]]
  p <- ggman(df,snp='cpg',bp='CpG_beg',chrom='CpG_chrm',pvalue='qval',sigLine=-log10(0.05), pointSize=2)
  # If any value passes FDR<0.05, color the point red
  if(any(df$qval<0.05)){
    ls <- df$cpg[df$qval<0.05]
    p <- ggmanHighlight(p,highlight=ls,colour = "red")
  }else{
    p <- p + scale_colour_grey(start = 0.5,end = 0.6)
  }
  # Add labels to top 5 CpGs (not used)
  # p <- ggmanLabel(p, labelDfm = head(df,5),snp='cpg',label='cpg',type = "text", colour="black")
  p <- p + coord_cartesian(ylim=c(0,6)) + ylab("-log10 FDR") + ggtitle(x) + theme_cowplot(10)
  return(p)
})
ggsave("plots/Manhattan_qval.pdf",width = 24,height = 12,
       plot_grid(plotlist = q_list[c(-1,-2)])
)

# Write filtered DMP list
lapply(names(dmp_man), function(x){
  dat <- dmp_man[[x]]
  colnames(dat)[4:5] <- c("chromosome_hg19","position_hg19")
  dat$chromosome_hg19 <- paste0("chr",dat$chromosome_hg19)
  write.table(dat,paste0("out/DMPs_filtered_",x,".tsv"),quote = F,sep = "\t",row.names = F)
})