
# attempt to analyze repertoire data like an RNA-seq experiment
# because we have presumably correlated individuals
# it looks like edgeR might not be sufficient for a proper treatment
# we're going to try to use a voom + limma combination

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
cont.matrix = makeContrasts(SPvsTHYinY = Y.SP - Y.THY,
                            SPvsTHYinO = O.SP - O.THY,
                            Diff = (O.SP - O.THY) - (Y.SP - Y.THY),
                            levels = design)

fit2 = contrasts.fit(fit,cont.matrix)
fit2 = eBayes(fit2)

# generate some lists of things
# young specific
y_tt = topTable(fit2,coef='SPvsTHYinY',number=Inf,adjust='fdr')

filename = paste('lists/','SPvsTHYinY_limma_top.txt',sep='')
write.table(y_tt,file=filename,row.names=FALSE,sep='\t')

# old specific
o_tt = topTable(fit2,coef='SPvsTHYinO',number=Inf,adjust='fdr')

filename = paste('lists/','SPvsTHYinO_limma_top.txt',sep='')
write.table(o_tt,file=filename,row.names=FALSE,sep='\t')

# difference
OvsY_diff_tt = topTable(fit2,coef='Diff',number=Inf,adjust='fdr')

filename = paste('lists/','OvsY_SPvsTHY_limma_top.txt',sep='')
write.table(OvsY_diff_tt,file=filename,row.names=FALSE,sep='\t')

# volcano some things
require(ggplot2)
volcano = function(tt,pre,special=NULL){
  g = ggplot(data=tt,mapping=aes(logFC,-log10(adj.P.Val)))
  g = g + geom_point(alpha=0.25,colour='blue')
  if (!is.null(special)){
    g = g + geom_point(mapping=aes(logFC,-log10(adj.P.Val)),data=tt[tt$Symbol %in% special,],colour='red')
    g = g + geom_text(mapping=aes(logFC,-log10(adj.P.Val),label=Symbol),data=tt[tt$Symbol %in% special,],hjust=1,vjust=1,colour='red')
  }
  g
  
  ggsave(paste('plots/',pre,'_volcano.png',sep=''))
}

volcano(y_tt,'y')
volcano(o_tt,'o')
volcano(OvsY_diff_tt,'OvsY')

# venn diagram to some general sense of whats going on 
results = decideTests(fit2,adjust.method='fdr',p.value=0.05)
vennDiagram(results)

lookup = function(rank,tt,IMGT.df){
  vf = IMGT.df[,'V_GENE'] == tt[rank,'V_GENE']
  jf = IMGT.df[,'J_GENE'] == tt[rank,'J_GENE']
  cf = IMGT.df[,'CDR3'] == tt[rank,'CDR3']
  return(IMGT.df[vf&jf&cf,])
}


