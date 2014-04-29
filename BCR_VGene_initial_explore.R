
# consolidate BCR data based on v gene families 

setwd("~/Documents/Dinner/Repertoire/3H9 BCR Sequencing")

require(plyr)

# load up full data
IMGT.df = read.table("IMGT Combined/IMGT_combined.txt",header=TRUE,stringsAsFactors=FALSE)
# only label by v gene
v_all_counts.df = IMGT.df[-(2:3)] # without labels

# aggregate based on the v gene column
# should end up with 123 different rows
vgene.df = ddply(v_all_counts.df, "V_GENE", numcolwise(sum))
vcounts.df = vgene.df[-1] # no label


# heatmap
# for a reasonable color scheme
library(gplots)
library(RColorBrewer)
bluescale = colorRampPalette(brewer.pal(9,"Blues"))(100)
# use count matrix to generate some heatmaps and things
cmat = sqrt(vcounts.df) # square root transform
ccor = cor(cmat)

pre = ''
png(paste('BCR V Gene plots/',pre,'all_present_heatmap.png',sep=''))
heatmap.2(ccor,Colv='Rowv',trace='none',density.info='none',col=bluescale)
dev.off()

