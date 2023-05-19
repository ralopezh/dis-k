import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

data = FortranFile('diagram.save','r')

nroots = data.read_ints(dtype=np.int32)[0]
nw = data.read_ints(dtype=np.int32)[0]
[k0, angle] = data.read_ints(dtype=float).reshape((2,nroots))
[wr, wi] = data.read_ints(dtype=float).reshape((2,nw))
zr = data.read_ints(dtype=float).reshape((nroots,nw,nw))
zi = data.read_ints(dtype=float).reshape((nroots,nw,nw))
iindx = data.read_ints(dtype=np.int32)[0]
cindx = chr(iindx) #convert int index into char

## Plot contour
wps = '\omega_{p'+cindx+'}'
Wcs = '\Omega_'+cindx

if nroots==1:
    fig, ax = plt.subplots()
    axs = [ax]
else:
    fig, axs = plt.subplots(1,nroots,figsize=(5*nroots,5),constrained_layout=True)

for i in range(0,nroots):
    axs[i].contour(wr, wi, zr[i,:,:], levels=[0.0], colors='blue')
    axs[i].contour(wr, wi, zi[i,:,:], levels=[0.0], colors='red')
    
    axs[i].set_xlabel(r'${\rm Re}(\omega/'+Wcs+')$')
    axs[i].set_ylabel(r'${\rm Im}(\omega/'+Wcs+')$')
    axs[i].set_title(r'$ck_0/'+wps+' = $'+str(k0[i])+r'$\,,\quad\theta = $'+str(angle[i]))

plt.show()

