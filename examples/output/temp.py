import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import matplotlib.colors as colors


data = FortranFile('temp.save','r')
[npts, nth, nroots] = data.read_ints(dtype=np.int32)
kvec = data.read_reals(dtype=float)
thvec = data.read_reals(dtype=float)
wsol = data.read_reals(dtype=complex).reshape((nroots,nth,npts))

wr = np.real(wsol)
wi = np.imag(wsol)

if nroots==1:
    fig, ax=plt.subplots(1, 1, figsize=(5,5),constrained_layout=True)
    axs = np.tile(ax,(1))
else:
    fig, axs = plt.subplots(nroots, 1, figsize=(5,5*nroots),constrained_layout=True)

nl = 35
for i in range(0,nroots):
    lmin2 = 1.0e-6 
    lmax2 = wi[i,:,:].max()
    dl2 = (lmax2-lmin2)/(nl-1)
    levs2 = np.arange(nl)*dl2+lmin2
    if lmin2<lmax2:
        cs2 = axs[i].contourf(kvec, thvec, wi[i,:,:], levels=levs2,cmap='jet')
        cbar = fig.colorbar(cs2,ax=axs[i],location='top',label=r'$\gamma$')
        cbar.ax.locator_params(nbins=5)
        axs[i].set_xlabel(r'$k$')
        axs[i].set_ylabel(r'$\theta$')
        
plt.show()