
# Poking around the repertoire data

setwd("~/Documents/Dinner/Repertoire/3H9 BCR Sequencing")
library(ggplot2)
library(reshape2)

IMGT.df = read.table("IMGT Combined/IMGT_combined.txt",header=TRUE,stringsAsFactors=FALSE)

counts.df = IMGT.df[-(1:3)]
combined_counts.df = cbind(apply(counts.df[1:5],1,sum),apply(counts.df[6:10],1,sum),
                                  apply(counts.df[11:15],1,sum),apply(counts.df[16:20],1,sum))
colnames(combined_counts.df) = c('Y.SP','Y.THY','O.SP','O.THY')

# get some general counts info
total_reads = apply(counts.df,2,sum)
unique_reads = apply(counts.df > 0,2,sum)
presence = apply(counts.df > 0,1,sum)
singleton_reads = apply(counts.df == 1,2,sum)
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

simple_bar(colnames(counts.df),total_reads,'plots/total_reads.png')
simple_bar(colnames(counts.df),unique_reads,'plots/unique_reads.png')
simple_bar(colnames(counts.df),singleton_reads,'plots/singleton_reads.png')
simple_bar(colnames(counts.df),P0s,'plots/P0s.png')

# separately deal with presence
zeros = numeric(20)
for (i in 1:20){
  zeros[i] = sum(presence==i)
}
simple_bar(1:20,zeros,'plots/presence.png')

# plot some rank plots
rank.df = data.frame()
for(j in 1:length(colnames(counts.df))){
  unique_vals = sort(unique(counts.df[counts.df[,j]>=1,j]),T)
  val_counts = lapply(unique_vals,function(x){sum(counts.df[,j] == x)})
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
  rank.df = rbind(rank.df,data.frame(sample=rep(colnames(counts.df)[j],length(unique_vals)),
                                     ranks=to_ranks,
                                     vals=to_vals))
}

# doesn't look like theres a particularly sharp decrease at any rank or threshold
g = ggplot(data=rank.df,aes(log10(ranks),log10(vals),group=sample,color=sample))
g = g + geom_line()
g

ggsave('plots/ranks.png')

# heatmap
# for a reasonable color scheme
library(gplots)
library(RColorBrewer)
bluescale = colorRampPalette(brewer.pal(9,"Blues"))(100)
# use count matrix to generate some heatmaps and things
cmat = sqrt(counts.df) # square root transform
ccor = cor(cmat)

pre = ''
png(paste('plots/',pre,'all_present_heatmap.png',sep=''))
heatmap.2(ccor,Colv='Rowv',trace='none',density.info='none',col=bluescale)
dev.off()

