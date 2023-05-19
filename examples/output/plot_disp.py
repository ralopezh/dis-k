import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

data = FortranFile('output.save','r')
[npts, nth, nsp, nroots, damped] = data.read_ints(dtype=np.int32)
[rpem, wpewce, wpwc] = data.read_ints(dtype=float)
[k0, th0] = data.read_reals(dtype=float).reshape((2,nroots))
[rm, ns, rq, betapal, anis, kappa, Ua] = data.read_ints(dtype=float).reshape((7,nsp))

kvec = data.read_reals(dtype=float)
wsol = data.read_reals(dtype=complex).reshape((nroots,npts))
pols = data.read_reals(dtype=complex).reshape((nroots,npts))
ratios = data.read_reals(dtype=float).reshape((8,nroots,npts))
iindx = data.read_ints(dtype=np.int32)[0]
cindx = chr(iindx) #convert int index into char

# extract real and imaginary part of omega
wr = np.real(wsol)
wi = np.imag(wsol)
pol = np.real(pols)

# Prepare ratios
rex2 = ratios[0,:,:]
rey2 = ratios[1,:,:]
rez2 = ratios[2,:,:]
rbx2 = ratios[3,:,:]
rby2 = ratios[4,:,:]
rbz2 = ratios[5,:,:]
R = ratios[6,:,:] #Ex/Ey
T = ratios[7,:,:] #Ez/Ey

## Plot result
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15,5),constrained_layout=True)
#plt.tight_layout()

for i in range(0,nroots):
    ax1.plot(kvec,wr[i])
    ax2.plot(kvec,wi[i])
    ax3.plot(kvec,pol[i])
#    ax3.plot(kvec,ratios[0,i]) # |Ex/E|^2
#    ax3.plot(kvec,ratios[1,i],linestyle='dashed') # |Ey/E|^2
#    ax3.plot(kvec,ratios[2,i],linestyle='dotted') # |Ez/E|^2
#    ax3.plot(kvec,ratios[3,i]) # |Bx/B|^2
#    ax3.plot(kvec,ratios[4,i],linestyle='dashed') # |By/B|^2
#    ax3.plot(kvec,ratios[5,i],linestyle='dotted') # |Bz/B|^2

wps = '\omega_{p'+cindx+'}'
Wcs = '\Omega_'+cindx

ax1.set_xlabel(r'$ck/'+wps+'$')
ax2.set_xlabel(r'$ck/'+wps+'$')
ax3.set_xlabel(r'$ck/'+wps+'$')
ax1.set_ylabel(r'$\omega/'+Wcs+'$')
ax2.set_ylabel(r'$\gamma/'+Wcs+'$')
ax3.set_ylabel(r'$P\,=\,{\Re}\{i{\rm Sign}(\omega_r)E_x/E_y\}$')
#
plt.show()
