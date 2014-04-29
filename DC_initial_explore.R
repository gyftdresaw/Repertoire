
# some initial looks at the DC TCR data

setwd("~/Documents/Dinner/Repertoire/3H9 BCR Sequencing")
library(ggplot2)
library(reshape2)

DC.df = read.table("TCR Combined/MBDC_Shared.txt",header=TRUE,stringsAsFactors=FALSE,sep='\t')
DC.df = DC.df[-18] # weird extra column showing up
DC.df = DC.df[c(1:3,14:17)]
present_in_DC = apply(DC.df[-(1:3)],1,sum)!=0

DC.df = DC.df[present_in_DC,]

DC_counts.df = DC.df[-(1:3)] # there's a lot, 430965 present in DC experiment

# get some general counts info
total_reads = apply(DC_counts.df,2,sum)
unique_reads = apply(DC_counts.df > 0,2,sum)
presence = apply(DC_counts.df > 0,1,sum)
singleton_reads = apply(DC_counts.df == 1,2,sum)
P0s = singleton_reads/total_reads # total probability of unseen species

# save some plots
simple_bar = function(cols,vals,to_file){
  # total reads
  g = ggplot(data=data.frame(x=cols,vals=vals),aes(x,vals))
  g = g + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle=90))
  # g = g + scale_y_log10()
  print(g)
  ggsave(to_file)
}

simple_bar(colnames(DC_counts.df),total_reads,'DC TCR plots/total_reads.png')
simple_bar(colnames(DC_counts.df),unique_reads,'DC TCR plots/unique_reads.png')
simple_bar(colnames(DC_counts.df),singleton_reads,'DC TCR plots/singleton_reads.png')
simple_bar(colnames(DC_counts.df),P0s,'DC TCR plots/P0s.png')

# separately deal with presence
zeros = numeric(4)
for (i in 1:4){
  zeros[i] = sum(presence==i)
}
simple_bar(1:4,zeros,'DC TCR plots/presence.png')

# plot some rank plots
rank.df = data.frame()
for(j in 1:length(colnames(DC_counts.df))){
  unique_vals = sort(unique(DC_counts.df[DC_counts.df[,j]>=1,j]),T)
  val_counts = lapply(unique_vals,function(x){sum(DC_counts.df[,j] == x)})
  val_counts = as.numeric(val_counts)
  
  to_ranks = numeric(2*length(unique_vals))
  to_vals = numeric(2*length(unique_vals))
  current_rank = 0
  for (i in 1:length(unique_vals)){
    current_rank = current_rank + 1 # increment to next rank
    to_ranks[2*i-1] = current_rank
    to_vals[2*i-1] = unique_vals[i]
    # update current_rank to end of region
    current_rank = current_rank + val_counts[i] - 1
    to_ranks[2*i] = current_rank
    to_vals[2*i] = unique_vals[i]
  }
  rank.df = rbind(rank.df,data.frame(sample=rep(colnames(DC_counts.df)[j],length(unique_vals)),
                                     ranks=to_ranks,
                                     vals=to_vals))
}

# doesn't look like theres a particularly sharp decrease at any rank or threshold
g = ggplot(data=rank.df,aes(log10(ranks),log10(vals),group=sample,color=sample))
g = g + geom_line()
g

ggsave('DC TCR plots/ranks.png')


# heatmap
# for a reasonable color scheme
library(gplots)
library(RColorBrewer)
bluescale = colorRampPalette(brewer.pal(9,"Blues"))(100)
# use count matrix to generate some heatmaps and things
cmat = sqrt(DC_counts.df) # square root transform
ccor = cor(cmat)

pre = ''
png(paste('DC TCR plots/',pre,'all_present_heatmap.png',sep=''))
heatmap.2(ccor,Colv='Rowv',trace='none',density.info='none',col=bluescale)
dev.off()







