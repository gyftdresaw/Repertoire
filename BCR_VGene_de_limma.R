
# differential expression on just the V Genes of the BCR experiment

require(edgeR)
require(limma)
require(plyr)

# load up full data
IMGT.df = read.table("IMGT Combined/IMGT_combined.txt",header=TRUE,stringsAsFactors=FALSE)
# only label by v gene
v_all_counts.df = IMGT.df[-(2:3)] # without labels

# aggregate based on the v gene column
# should end up with 123 different rows
vgene.df = ddply(v_all_counts.df, "V_GENE", numcolwise(sum))
vcounts.df = vgene.df[-1] # no label

# labeling the sample groups/types
ages = c(rep('Y',10),rep('O',10))
types = c(rep('SP',5),rep('THY',5),rep('SP',5),rep('THY',5))
individuals = c(rep(1:5,2),rep(6:10,2)) # unique labels for individuals
conditions = paste(ages,types,sep='.') # combined unique conditions

# data normalization 
# first construct a DGEList to contain the count/annotation information
vgene_info = vgene.df[1]

dge_data = DGEList(counts=vcounts.df,genes=vgene_info)
# normalize this dge data
# method is 'TMM', supposedly accounting for
# composition differences in RNAseq (overrepresentation of some samples?)
# seems to assume that a majority of genes are not differentially expressed
# which is not necessarily true in this context
# (maybe should just use simple good turing estimation instead of normalization?)
dge_data = calcNormFactors(dge_data)

# create appropriate design matrix
# we're going to try to keep all the coefficients as explicit as possible
design = model.matrix(~0+factor(conditions))
colnames(design) = levels(factor(conditions))

# transform to logCPM (log counts per million)
# should be log2, calculated with offset of 0.5 
# which is arbitrary but probably something reasonable to do
logCPM_data = voom(dge_data,design,plot=TRUE)
# whether or not the data was normalized as above doesn't 
# really change the mean-variance trend, which doesn't look so great
# the trend should really decrease and then level off at higher counts
# but seems to increase artificially because of the low counts and then
# decrease slowly at higher counts.
#
# this looks slightly better behave with presence > 5

# pairs on logCPM_data doesn't look too different from pairs
# directly on log'd count data or the sgt proportions

# now we move the analysis to limma
# i think this is taking an unreasonably long time with every candidate
# going to filter somewhere above
# with presence > 5 still kinda slow for some reason
dupcor = duplicateCorrelation(logCPM_data,design,block=factor(individuals))
# dupcor$consensus pretty low (0.11)

# fit linear model using this individual duplicate correlation
fit = lmFit(logCPM_data, design, block=factor(individuals), correlation=dupcor$consensus)

# construct contrast matrix for comparisons
cont.matrix = makeContrasts(SPvsTHYinY = Y.SP - Y.THY,
                            SPvsTHYinO = O.SP - O.THY,
                            Diff = (O.SP - O.THY) - (Y.SP - Y.THY),
                            levels = design)

fit2 = contrasts.fit(fit,cont.matrix)
fit2 = eBayes(fit2)

# generate some lists of things
# young specific
vy_tt = topTable(fit2,coef='SPvsTHYinY',number=Inf,adjust='fdr')

filename = paste('BCR V Gene lists/','SPvsTHYinY_vgene_limma_top.txt',sep='')
write.table(vy_tt,file=filename,row.names=FALSE,sep='\t')

# old specific
vo_tt = topTable(fit2,coef='SPvsTHYinO',number=Inf,adjust='fdr')

filename = paste('BCR V Gene lists/','SPvsTHYinO_vgene_limma_top.txt',sep='')
write.table(vo_tt,file=filename,row.names=FALSE,sep='\t')

# difference
vOvsY_diff_tt = topTable(fit2,coef='Diff',number=Inf,adjust='fdr')

filename = paste('BCR V Gene lists/','OvsY_SPvsTHY_vgene_limma_top.txt',sep='')
write.table(vOvsY_diff_tt,file=filename,row.names=FALSE,sep='\t')

######################
# heatmap expression #
######################

# heatmap
# for a reasonable color scheme
library(gplots)
library(RColorBrewer)
bluescale = colorRampPalette(brewer.pal(9,"Blues"))(100)

# create big data matrix of log count data
log_vmat = logCPM_data$E
rownames(log_vmat) = vgene.df[,1]

pre = ''
png(paste('BCR V Gene plots/',pre,'v_genes_heatmap.png',sep=''),width=1000,height=2000)
heatmap.2(log_vmat,Colv=FALSE,trace='none',density.info='none',scale='none',
          col=redgreen(75),srtCol=90,cexRow=1.5,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()

# seems to be mostly grouping by net expression 
# subtract each row by the mean? 
row_scaled_vmat = t(scale(t(log_vmat)))

pre = ''
png(paste('BCR V Gene plots/',pre,'v_genes_row_scaled_heatmap.png',sep=''),width=1000,height=2000)
heatmap.2(row_scaled_vmat,Colv=FALSE,trace='none',density.info='none',scale='none',
          col=redgreen(75),srtCol=90,cexRow=1.5,cexCol=1.5,margins=c(9,9),keysize=0.5)
dev.off()








