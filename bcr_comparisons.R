
# basically perform de_limma but for some specific BCR comparisons
# 
# comparisons are:
# YSPxYTHY, OSPxOTHY, YTHYxOTHY, and YSPxOSP
#
# also make some nicer lists

# starting off with same load data and conversion as de_limma

require(edgeR)
require(limma)

setwd("~/Documents/Dinner/Repertoire/3H9 BCR Sequencing")

# load em up

IMGT.df = read.table("IMGT Combined/IMGT_combined.txt",header=TRUE,stringsAsFactors=FALSE)
counts.df = IMGT.df[-(1:3)] # without labels

# labeling the sample groups/types
ages = c(rep('Y',10),rep('O',10))
types = c(rep('SP',5),rep('THY',5),rep('SP',5),rep('THY',5))
individuals = c(rep(1:5,2),rep(6:10,2)) # unique labels for individuals
conditions = paste(ages,types,sep='.') # combined unique conditions

# data normalization 
# first construct a DGEList to contain the count/annotation information
gene_info = IMGT.df[1:3]

# we're going to filter right here
# only keep things that have reads in more than 5 samples
presence = apply(counts.df > 0,1,sum)

# new filtered data
# having at least one read in more than 5 samples trims the list to 35707 combinations
fcounts.df = counts.df[presence > 5,]
fgene_info = gene_info[presence > 5,]


dge_data = DGEList(counts=fcounts.df,genes=fgene_info)
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
# YSPxYTHY, OSPxOTHY, YTHYxOTHY, and YSPxOSP
cont.matrix = makeContrasts(THYvsSPinY = Y.THY - Y.SP,
                            THYvsSPinO = O.THY - O.SP,
                            OTHYvsYTHY = O.THY - Y.THY,
                            OSPvsYSP = O.SP - Y.SP,
                            levels = design)

fit2 = contrasts.fit(fit,cont.matrix)
fit2 = eBayes(fit2)

### get the top table results ###
TvSinY_tt = topTable(fit2,coef='THYvsSPinY',number=Inf,adjust='fdr')
TvSinO_tt = topTable(fit2,coef='THYvsSPinO',number=Inf,adjust='fdr')
OTvYT_tt = topTable(fit2,coef='OTHYvsYTHY',number=Inf,adjust='fdr')
OSPvYSP_tt = topTable(fit2,coef='OSPvsYSP',number=Inf,adjust='fdr')

# make some nice lists
# to include:
# V,J,CDDR3, fold change, ave reads in groups, pval

# need to get ave reads in each group
ave_fcounts.df = data.frame(Y.SP = rowMeans(fcounts.df[,conditions=='Y.SP']),
                            Y.THY = rowMeans(fcounts.df[,conditions=='Y.THY']),
                            O.SP = rowMeans(fcounts.df[,conditions=='O.SP']),
                            O.THY = rowMeans(fcounts.df[,conditions=='O.THY']))
# these counts correspond to V,J,CDR3 in fgene_info
# need identifier for these
fgene_id = apply(fgene_info,1,function(x){paste(x,collapse='.')})

# create the lists

# YSPxYTHY
dtt = TvSinY_tt
dtt_id = apply(dtt[,1:3],1,function(x){paste(x,collapse='.')})
TvSinY_list = data.frame(V_GENE=dtt$V_GENE,J_GENE=dtt$J_GENE,CDR3=dtt$CDR3,
                         FOLD_CHANGE=sign(dtt$logFC) * 2^abs(dtt$logFC),
                         AVE_YTHY=ave_fcounts.df[match(dtt_id,fgene_id),'Y.THY'],
                         AVE_YSP=ave_fcounts.df[match(dtt_id,fgene_id),'Y.SP'],
                         FDR_PVAL=dtt$adj.P.Val,UNADJUSTED_PVAL=dtt$P.Value)

# OSPxOTHY
dtt = TvSinO_tt
dtt_id = apply(dtt[,1:3],1,function(x){paste(x,collapse='.')})
TvSinO_list = data.frame(V_GENE=dtt$V_GENE,J_GENE=dtt$J_GENE,CDR3=dtt$CDR3,
                         FOLD_CHANGE=sign(dtt$logFC) * 2^abs(dtt$logFC),
                         AVE_OTHY=ave_fcounts.df[match(dtt_id,fgene_id),'O.THY'],
                         AVE_OSP=ave_fcounts.df[match(dtt_id,fgene_id),'O.SP'],
                         FDR_PVAL=dtt$adj.P.Val,UNADJUSTED_PVAL=dtt$P.Value)

# YTHYxOTHY
dtt = OTvYT_tt
dtt_id = apply(dtt[,1:3],1,function(x){paste(x,collapse='.')})
OTvYT_list = data.frame(V_GENE=dtt$V_GENE,J_GENE=dtt$J_GENE,CDR3=dtt$CDR3,
                        FOLD_CHANGE=sign(dtt$logFC) * 2^abs(dtt$logFC),
                        AVE_OTHY=ave_fcounts.df[match(dtt_id,fgene_id),'O.THY'],
                        AVE_YTHY=ave_fcounts.df[match(dtt_id,fgene_id),'Y.THY'],
                        FDR_PVAL=dtt$adj.P.Val,UNADJUSTED_PVAL=dtt$P.Value)

# YSPxOSP
dtt = OSPvYSP_tt
dtt_id = apply(dtt[,1:3],1,function(x){paste(x,collapse='.')})
OSPvYSP_list = data.frame(V_GENE=dtt$V_GENE,J_GENE=dtt$J_GENE,CDR3=dtt$CDR3,
                         FOLD_CHANGE=sign(dtt$logFC) * 2^abs(dtt$logFC),
                         AVE_OSP=ave_fcounts.df[match(dtt_id,fgene_id),'O.SP'],
                         AVE_YSP=ave_fcounts.df[match(dtt_id,fgene_id),'Y.SP'],
                         FDR_PVAL=dtt$adj.P.Val,UNADJUSTED_PVAL=dtt$P.Value)

# write these lists to file
# TvSinY_list, TvSinO_list, OTvYT_list, OSPvYSP_list
filename = paste('BCR comparison lists/','YTHYvsYSP_limma.txt',sep='')
write.table(TvSinY_list,file=filename,row.names=FALSE,sep='\t')

filename = paste('BCR comparison lists/','OTHYvsOSP_limma.txt',sep='')
write.table(TvSinO_list,file=filename,row.names=FALSE,sep='\t')

filename = paste('BCR comparison lists/','OTHYvsYTHY_limma.txt',sep='')
write.table(OTvYT_list,file=filename,row.names=FALSE,sep='\t')

filename = paste('BCR comparison lists/','OSPvsYSP_limma.txt',sep='')
write.table(OSPvYSP_list,file=filename,row.names=FALSE,sep='\t')

# volcano some things
require(ggplot2)
volcano = function(tt,pre){
  g = ggplot(data=tt,mapping=aes(logFC,-log10(adj.P.Val)))
  g = g + geom_point(alpha=0.25,colour='blue')
  g = g + xlab('Log2(Fold Change)') + ylab('-Log10(FDR P Value)')
  g
  
  ggsave(paste('BCR comparison plots/',pre,'_volcano.png',sep=''))
}

volcano(TvSinY_tt,'YTHYvsYSP')
volcano(TvSinO_tt,'OTHYvsOSP')
volcano(OTvYT_tt,'OTHYvsYTHY')
volcano(OSPvYSP_tt,'OSPvsYSP')



