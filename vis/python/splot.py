import sys
import numpy as np
import matplotlib.pyplot as plt

dir      = './'
filename ='species_g0_r1_z0.out'

name     = []
weight   = []
with open(dir+filename) as hf:
  for line in hf:
    data   = line.split()
    weight = np.append(weight,float(data[0]))
    name   = np.append(name,data[1])

x = np.linspace(1,len(weight),len(weight))

plt.scatter(x,weight,'k',linewidth=2)
plt.axhline(y=1.0)
plt.xlabel('species')
plt.ylabel('weight')

plt.axis([0,len(weight),1e-5,1e1])
plt.show()
