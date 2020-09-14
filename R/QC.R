options(stringsAsFactors = F)

###################
## Load raw data ##
###################

library(minfi)
targets <- read.metharray.sheet("Data")
RGSet <- read.metharray.exp(targets=targets)

########
## QC ##
########

## Generate qc report
qcReport(RGSet,sampGroups=targets$Slide,pdf="qcReport.pdf")

## Generate intensity plot
Mset <- preprocessRaw(RGSet)
pdf("plots/QC_450k_intensity.pdf")
plotQC(getQC(Mset))
dev.off()

## Plot bisulfite conversion
library(shinyMethyl)
RGSet.SM <- shinySummarize(RGSet)
#runShinyMethyl(RGSet.SM)

# Function adapted from ShinyMethyl source code
bisulfite.plot <- function(x, ... ){
  require(shinyMethyl)
  require(ggplot2)
  require(cowplot)
  grn <- getGreenControls(x)
  red <- getRedControls(x)
  
  bc1 <- (apply(grn[["BISULFITE CONVERSION I"]][1:3,], 2, mean) + apply(red[["BISULFITE CONVERSION I"]][7:9,], 2, mean))/2 
  bc2 <- apply(red[["BISULFITE CONVERSION II"]], 2, mean)
  
  gridExtra::grid.arrange(nrow=1,
    qplot(y=bc1, x=sampleNames(RGSet.SM), ... ) + coord_cartesian(ylim=c(0,4e4)) + labs(title="Bisulfite Conversion Probe type I", y="Control intensity") + theme_cowplot(10) + theme(axis.title.y = element_blank()) + coord_flip() + guides(col=F),
    qplot(y=bc2, x=sampleNames(RGSet.SM), ... ) + coord_cartesian(ylim=c(0,4e4)) + labs(title="Bisulfite Conversion Probe type II", y="Control intensity") + theme_cowplot(10) + theme(axis.title.y = element_blank()) + coord_flip() + guides(col=F)
  )
}

library(ggplot2)
library(cowplot)
ggsave("plots/QC_bisulfiteConversion.pdf",width = 12,height = 8,
  bisulfite.plot(RGSet.SM, col=targets$Slide)
)

## Detection P-values
detP <- detectionP(RGSet)
failed <- detP>0.01

ggsave("plots/QC_detectionP.pdf",
  ggplot(as.data.table(melt(colMeans(failed)),keep.rownames = T), aes(x=rn, y=value)) + geom_point() + labs(x=NULL, y= "Fraction of failed positions") + coord_flip() + theme_cowplot(10)
)

## Plot density
# Generates ggplot2 color palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

pdf("plots/QC_density.pdf")
densityPlot(RGSet,sampGroups = targets$Slide,pal=gg_color_hue(4))
lines(density(getBeta(RGSet)[,"100896160012_R04C01"],na.rm=T),lwd=3,col=gg_color_hue(1))
text(0.4,1,"100896160012_R04C01", col=gg_color_hue(1))
dev.off()