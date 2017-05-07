import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 16}
mpl.rc('font', **font) 

height = 4.
col = ['r','k','b','g','c','m',\
      'y','#8A2BE2','#DEB888','#FF8C00','#808080']
# read data
r = [1,10,100]
for i in range(0,3):
  print (r[i],'AU')
  with open("scplx"+str(r[i]), "r") as file:
    res = [[float(x) for x in line.split()] for line in file]
  cplx = np.array(res)
  with open("ssmpl"+str(r[i]), "r") as file:
    res = [[float(x) for x in line.split()] for line in file]
  smpl = np.array(res)

  z = np.linspace(0,height,cplx.shape[0])
  z = np.logspace(-1,-4,cplx.shape[0])
  for j in range(cplx.shape[1]-1):
    ind = cplx.shape[1]-1
    plt.semilogx(z,np.log10(cplx[:,j+1]/cplx[:,0]),linestyle='-',linewidth=1.5,color=col[ind*i+j],label=r'{}AU'.format(r[i]))
    plt.semilogx(z,np.log10(smpl[:,j+1]/smpl[:,0]),linestyle='--',linewidth=1.5,color=col[ind*i+j])

plt.legend(loc=4)
plt.xlabel(r'$\Sigma[g\cdot cm^2]$')
plt.ylabel(r'$\log(n/n_H)$')
plt.axis([z[0],z[-1],-7,-2])
#fig = mpl.pyplot.gcf()
#fig.set_size_inches(18,6)
plt.savefig('senum.pdf',format='pdf',dpi=300,bbox_inches='tight')
plt.show()
