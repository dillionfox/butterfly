import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

import sys

# This program takes a tab delimited file of the form "NN    #####" where the first column contains
# an alphabetical list of amino acid subsequences of length 2 and the second column contains an
# integer value counting the number of occurences of that subsequence in a given dataset and plots this 
# on a heatmap

xax = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
yax = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

idx = range(0,19)

trans = zip(xax, idx) 

count = np.zeros((len(xax), len(yax)))

for line in open(sys.argv[1],'r'):
    for data in line.splitlines():
        x, y = 0, 0 
        aa, ct = data.strip().split()
        for elem in trans:
            if elem[0] == aa[0]:
                x = elem[1]
        for elem in trans:
            if elem[0] == aa[1]:
                y = elem[1]
        count[x,y] = ct

fig, ax = plt.subplots()
#im = ax.imshow(count, interpolation='none')

sns.heatmap(count, cmap='Greys_r',square=True, linewidths=.5)

ax.set_xticklabels(xax)
ax.set_yticklabels(yax)

ax.xaxis.tick_top()


plt.yticks(rotation=0)

plt.tick_params(which='both', bottom=False, top=False, left=False, right=False, labelbottom=False)
ax.tick_params(axis='x', which='major', pad=0)

fig.set_tight_layout(True)

plt.savefig(sys.argv[1]+'.png')
