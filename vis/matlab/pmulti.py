import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 16}
mpl.rc('font', **font) 
fig,ax = plt.subplots(1,3,sharey= True,sharex=True)

height = 4.
sp = ['H2O','CO','NH3','H2']
sp = ['HCO+','NH4+','H3+','H3O+']
sp = ['C+','C','CO']
sp = ['S','S+','SO','SO+']
col = ['r','k','b','g','c','m',\
      'y','#8A2BE2','#DEB888','#FF8C00','#808080']
#sp = ['e-','gr[2-]','gr[-]','gr','gr[+]']

# read data
r = [1,10,100]
for i in range(0,3):
  with open("scplxm"+str(r[i]), "r") as file:
    res = [[float(x) for x in line.split()] for line in file]
  cplx = np.array(res)
  with open("ssmplm"+str(r[i]), "r") as file:
    res = [[float(x) for x in line.split()] for line in file]
  smpl = np.array(res)

  #z = np.linspace(0,height,cplx.shape[0])
  if i==0:
    z = np.logspace(-1,-4,cplx.shape[0])
  for j in range(cplx.shape[1]-1):
    ax[i].semilogx(z,np.log10(1e-50+cplx[:,j+1]/cplx[:,0]),linestyle='-',linewidth=1.5,color=col[j],label=sp[j])
    print (z[0])
    ax[i].semilogx(z,np.log10(1e-50+smpl[:,j+1]/smpl[:,0]),linestyle='--',linewidth=1.5,color=col[j])
    ax[i].set_xlabel(r'$\Sigma[g\cdot cm^2]$')
    ax[i].set_title(r'{}AU'.format(r[i]))

ax[2].legend(loc=4,prop={'size':12})
ax[0].set_ylabel(r'$\log (n/n_H)$')
plt.axis([z[0],z[-1],-14,-4])
fig = mpl.pyplot.gcf()
fig.set_size_inches(18,4)
plt.savefig('ssden.pdf',format='pdf',dpi=300,bbox_inches='tight')
plt.show()
