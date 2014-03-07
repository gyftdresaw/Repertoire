
# do limma differential expression analysis on TCR data

require(edgeR)
require(limma)

setwd("~/Documents/Dinner/Repertoire/3H9 BCR Sequencing")

# read it up
TCR.df = read.table("TCR Combined/MB_Shared.txt",header=TRUE,stringsAsFactors=FALSE,sep='\t')
TCR.df = TCR.df[-14] # weird extra column showing up

TCR_counts.df = TCR.df[-(1:3)]

# labeling the sample groups/types
types = c(rep('WT',5),rep('KO',5))

# data normalization 
# first construct a DGEList to contain the count/annotation information
gene_info = TCR.df[1:3]

# we're going to filter right here
# only keep things that have reads in more than 5 samples
presence = apply(TCR_counts.df > 0,1,sum)

# new filtered data, down to 21823 combinations
fcounts.df = TCR_counts.df[presence > 5,]
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
# we only have two cases
design = model.matrix(~0+factor(types))
colnames(design) = levels(factor(types))

# transform to logCPM (log counts per million)
# should be log2, calculated with offset of 0.5 
# which is arbitrary but probably something reasonable to do
logCPM_data = voom(dge_data,design,plot=TRUE)
# whether or not the data was normalized as above doesn't 
# really change the mean-variance trend, which doesn't look so great
# the trend should really decrease and then level off at higher counts
# but seems to increase artificially because of the low counts and then
# decrease slowly at higher counts.

# limma without individual correlation
fit = lmFit(logCPM_data, design)

# construct contrast matrix for comparisons
cont.matrix = makeContrasts(KOvsWT = KO - WT,
                            levels = design)

fit2 = contrasts.fit(fit,cont.matrix)
fit2 = eBayes(fit2)

# generate some lists of things
# young specific
tcr_tt = topTable(fit2,coef='KOvsWT',number=Inf,adjust='fdr')

filename = paste('TCR lists/','KOvsWT_limma_top.txt',sep='')
write.table(tcr_tt,file=filename,row.names=FALSE,sep='\t')

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
  
  ggsave(paste('TCR plots/',pre,'_volcano.png',sep=''))
}

volcano(tcr_tt,'tcr')

# venn diagram to some general sense of whats going on 
results = decideTests(fit2,adjust.method='fdr',p.value=0.20)
vennDiagram(results)

lookup = function(rank,tt,IMGT.df){
  vf = IMGT.df[,'V.gene'] == tt[rank,'V.gene']
  jf = IMGT.df[,'J.gene'] == tt[rank,'J.gene']
  cf = IMGT.df[,'CDR3'] == tt[rank,'CDR3']
  return(IMGT.df[vf&jf&cf,])
}



