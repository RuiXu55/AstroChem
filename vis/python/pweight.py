import sys
import numpy as np
import matplotlib.pyplot as plt
import module

# construct a new Speciesname output
weight = []
with open('allspecies.out') as f:
  for line in f:
    data = line.split()
    weight = np.append(weight,np.float(data[0]))

ind = np.linspace(1,len(weight),len(weight))
plt.scatter(ind, weight,marker='d', s=28, c='b', alpha=0.6)
plt.xlabel('Species Index')
plt.ylabel('B[i]')
plt.yscale('log')
#plt.axhline(0.4,color='k',linestyle='--',linewidth=2)
plt.axis([0,100,1e-6,1e1])
plt.savefig('weight.pdf',format='pdf',dpi=300,bbox_inches='tight')
plt.show()
