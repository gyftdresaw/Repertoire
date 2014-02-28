
# try out some Simple Good Turing frequency estimate
# based on Gale and Sampson's (1995/2001) "Simple 
# Good Turing" algorithm

# implementation done by Max Bane
# github.com/maxbane/simplegoodturing

'''
REFERENCES:
    William Gale and Geoffrey Sampson. 1995. Good-Turing frequency estimation
    without tears. Journal of Quantitative Linguistics, vol. 2, pp. 217--37.
    
    See also the corrected reprint of same on Sampson's web site.
'''

import numpy as np
import pandas as pd

# load up consolidated repertoire data
# we're going to try pandas for fun
# should probably be like R
rep_data = pd.read_table('IMGT Combined/IMGT_combined.txt')

'''
# iterate through data frame and construct dictionary of counts
Nrows = rep_data.shape[0]
Nsamples = 20
count_dicts = [{} for i in xrange(Nsamples)] # list of empty dictionaries

# this is fairly slow
for i in xrange(Nrows):
    for j in xrange(1):
        # first three columns of data frame are ID
        # we won't actually keep track of the ID tuple
        # but instead just index by an int
        count_dicts[j][i] = rep_data.iloc[i,3+j]

# trying to convert one of these dictionaries into 
# a 'count of counts' dictionary via the sgt code
# is really slow, should probably try to get by without
# these dictionaries

# maybe should try this on the shakespeare case
# to test it out first
'''

# let's try to directly construct the count of counts dictionary
# for the sample indexed by i, retrieve the count of counts dictionary
# as well as the normal counts dictionary
def get_count_dicts(i):
    # get list of possible values
    vals = set(rep_data.iloc[:,3+i])
    if 0 in vals:
        vals.remove(0) # neglect 0 counts

    print 'got vals'
    # turn set into sorted list
    sorted_vals = sorted(list(vals))
    # for each one of these vals, need to count how many instances there were
    # takes a few moments but is manageable
    nonzero_data = rep_data[rep_data.iloc[:,3+i] != 0]
    print 'counting vals'
    counted_vals = [sum(nonzero_data.iloc[:,3+i]==sorted_vals[j]) for j in xrange(len(sorted_vals))]
    print 'done counting vals'
    count_dict = dict(zip(sorted_vals,counted_vals)) # hold count of counts
    # we'll also want to provide the normal count dictionary mapping species to counts
    species = [(nonzero_data.iloc[j,0],nonzero_data.iloc[j,1],nonzero_data.iloc[j,2]) 
               for j in xrange(len(nonzero_data))]
    counts_of_species = list(nonzero_data.iloc[:,3+i])
    counts = dict(zip(species,counts_of_species))

    return counts,count_dict

# try simple good turing estimation
# sometimes useful for autoreloading modules
# %load_ext autoreload
# %autoreload 2
import simplegoodturingmaster.sgt_HGMOD as sgt
import matplotlib.pyplot as plt
import numpy as np

def sgt_analysis(counts,count_dict,i):
    # do the sgt estimation
    sgtProbs,p0,reg_tuple = sgt.simpleGoodTuringProbs(counts,count_dict)
    
    # plot the results
    # look at comparison plot
    plt.figure()
    sgt.plotFreqVsGoodTuring(counts,count_dict,1.96,True)
    plt.savefig('plots/compare_freqs_%d.png' % i)

    rs,zs,a,b = reg_tuple
    # plot log10(count of counts) vs log10(counts) to see appropriate
    # log linear relationship
    plt.figure()

    # first plot real count data
    sorted_vals,counted_vals = zip(*count_dict.items())
    ax1 = plt.subplot(121)
    plt.scatter(sorted_vals,counted_vals,c='blue',alpha=0.2)
    ax1.set_yscale('log', basey=10)
    ax1.set_xscale('log', basex=10)
    ax1.set_xlabel('log10(r),frequency')
    ax1.set_ylabel('log10(Nr),frequency of frequency')

    # plot interpolation and linear fit
    ax2 = plt.subplot(122,sharex=ax1,sharey=ax1)
    plt.scatter(rs,zs,c='blue',alpha=0.2)
    plt.plot(rs,np.exp(a*np.log(np.array(rs)) + b),c='red')
    ax2.set_yscale('log', basey=10)
    ax2.set_xscale('log', basex=10)
    ax2.set_xlabel('log10(r),frequency')
    ax2.set_ylabel('log10(Zr),frequency of frequency')

    plt.suptitle('Regression: log(z) = %.2f*log(r) + %.2f' % (a,b))
    plt.tight_layout()
    plt.savefig('plots/interpolation_fit_%d.png' % i)

    # return results from sgt estimation
    return sgtProbs, p0, reg_tuple

# for each sample, we want to:
# retrieve count dictionaries
# perform sgt estimation
# reconstruct data based on sgt estimation

# perform data reconstruction
# get full list of species 
# kind of slow, only do this once
full_species = [(rep_data.iloc[j,0],rep_data.iloc[j,1],rep_data.iloc[j,2]) 
           for j in xrange(len(rep_data))]

sorted_full_species = sorted(full_species)

# dictionary to map species to appropriate index
species_indices = dict(zip(sorted_full_species,range(len(sorted_full_species))))

# take sgtProbs dictionary mapping species to probs 
# and return full probability vector based on sorted full species array
def get_full_probs(sgtProbs):
    species_probs = np.zeros(len(sorted_full_species))
    for k,v in sgtProbs.iteritems():
        species_probs[species_indices[k]] = v
    return species_probs

NSamples = 20
prob_data = np.empty((len(sorted_full_species),NSamples))
for i in xrange(NSamples):
    # get count dicts and perform analysis
    counts,count_dict = get_count_dicts(i)
    sgtProbs, p0, reg_tuple = sgt_analysis(counts,count_dict,i)

    prob_data[:,i] = get_full_probs(sgtProbs)

# write prob data to file
import csv

to_file = open('IMGT Probs/IMGT_probs.txt', 'w')
writer = csv.writer(to_file,delimiter='\t')

# write header
column_names = list(rep_data.columns)
writer.writerow(column_names[:3+NSamples])
for i in xrange(len(sorted_full_species)):
    writer.writerow(list(sorted_full_species[i]) + list(prob_data[i,:]))

to_file.close()


