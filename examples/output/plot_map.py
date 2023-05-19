import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import matplotlib.colors as colors

data = FortranFile('output_kth.save','r')
[npts, nth, nsp, nroots, damped] = data.read_ints(dtype=np.int32)
[rpem, wpewce, wpwc] = data.read_reals(dtype=float)
[k0, th0] = data.read_reals(dtype=float).reshape((2,nroots))
[rm, ns, rq, betapal, anis, kappa, Ua] = data.read_reals(dtype=float).reshape((7,nsp))

kvec = data.read_reals(dtype=float)
thvec = data.read_reals(dtype=float)
wsol = data.read_reals(dtype=complex).reshape((nroots,nth,npts))
pols = data.read_reals(dtype=complex).reshape((nroots,nth,npts))
ratios = data.read_reals(dtype=float).reshape((8,nroots,nth,npts))
iindx = data.read_ints(dtype=np.int32)[0]
cindx = chr(iindx) #convert int index into char

# extract real and imaginary part of omega
wr = np.real(wsol)
wi = np.imag(wsol)
pol = np.real(pols)

## Plot contour
if nroots==1:
    fig, ax=plt.subplots(1, 3, figsize=(15,5),constrained_layout=True)
    axs = np.tile(ax,(1,1))
else:
    fig, axs = plt.subplots(nroots, 3, figsize=(15,5*nroots),constrained_layout=True)

#levels
nl = 35

#labels
wps = '\omega_{p'+cindx+'}'
Wcs = '\Omega_'+cindx

#contours
for i in range(0,nroots):
    lmin1 = wr[i,:,:].min()
    lmax1 = wr[i,:,:].max()
    dl1 = (lmax1-lmin1)/(nl-1)
    levs1 = np.arange(nl)*dl1+lmin1
    if lmin1<lmax1:
        cs1 = axs[i,0].contourf(kvec, thvec, wr[i,:,:], levels=levs1,cmap='jet')
        cbar = fig.colorbar(cs1,ax=axs[i,0],location='top',label=r'$\omega_r/'+Wcs+'$')
        cbar.ax.locator_params(nbins=5)
        axs[i,0].set_xlabel(r'$ck/'+wps+'$')
        axs[i,0].set_ylabel(r'$\theta$')
    
for i in range(0,nroots):
    lmin2 = 1.0e-4 if damped==0 else wi[i,:,:].min()
    lmax2 = wi[i,:,:].max()
    dl2 = (lmax2-lmin2)/(nl-1)
    levs2 = np.arange(nl)*dl2+lmin2
    if lmin2<lmax2:
        cs2 = axs[i,1].contourf(kvec, thvec, wi[i,:,:], levels=levs2,cmap='jet')
        cbar = fig.colorbar(cs2,ax=axs[i,1],location='top',label=r'$\gamma/'+Wcs+'$')
        cbar.ax.locator_params(nbins=5)
        axs[i,1].set_xlabel(r'$ck/'+wps+'$')
        axs[i,1].set_ylabel(r'$\theta$')
        
for i in range(0,nroots):
    lmin3 = pol[i,:,:].min()
    lmax3 = pol[i,:,:].max()
    dl3 = (lmax3-lmin3)/(nl-1)
    levs3 = np.arange(nl)*dl3+lmin3
    if lmin3<lmax3:
        cs3 = axs[i,2].contourf(kvec, thvec, pol[i,:,:], levels=levs3,cmap='jet')
        cbar = fig.colorbar(cs3,ax=axs[i,2],location='top',\
        label=r'$P\,=\,{\rm Re}\{i{\rm Sign}(\omega_r)E_x/E_y\}$')
        cbar.ax.locator_params(nbins=5)
        axs[i,2].set_xlabel(r'$ck/'+wps+'$')
        axs[i,2].set_ylabel(r'$\theta$')
        

plt.show()
