"""
	Check gmapi indices for interpolation GLBb0.08 --> IAS 0.03
"""
from IPython import get_ipython
get_ipython().magic('reset -sf')

import os
import numpy as np
from copy import copy
import importlib
import matplotlib.pyplot as plt
import sys

sys.path.append('/home/ddmitry/codes/MyPython/hycom_utils')
sys.path.append('/home/ddmitry/codes/MyPython/draw_map')
sys.path.append('/home/ddmitry/codes/MyPython')

from mod_utils_fig import bottom_text

pthias = '/Net/kronos/ddmitry/hycom/TSIS/topo_IAS0.03/'
ftopo_ias = 'regional.depth'
fgrid_ias = 'regional.grid'
ftopo_sub = 'depth_09m11_ias0.03'

import mod_read_hycom
importlib.reload(mod_read_hycom)
from mod_read_hycom import read_grid_topo

LON, LAT, HH = read_grid_topo(pthias,ftopo_ias,fgrid_ias)
HH[HH<1.e10] = -1.*HH[HH<1.e10]
mm = HH.shape[0]
nn = HH.shape[1]

LON,LAT,HHglb = read_grid_topo(pthias,ftopo_sub,fgrid_ias)
HHglb[HHglb<1.e10] = -1.*HHglb[HHglb<1.e10]

clrmp = plt.cm.rainbow
clrmp.set_over(color=[0.4,0.4,0.4])

plt.ion()  # enables interactive mode
fig1 = plt.figure(1,figsize=(8,8))
plt.clf()
ax1 = plt.subplot(2,1,1)
im1 = plt.pcolormesh(LON,LAT,HH,shading='flat',cmap=clrmp,vmin=-8000., vmax=0)
ax1.axis('scaled')
ax1.set_xlim([np.min(LON), np.max(LON)])
ax1.set_ylim([np.min(LAT), np.max(LAT)])
plt.colorbar(im1)
plt.title('HH IAS '+ftopo_ias)

ax2 = plt.subplot(2,1,2)
im2 = plt.pcolormesh(LON,LAT,HHglb,shading='flat',cmap=clrmp,vmin=-8000., vmax=0)
ax2.axis('scaled')
ax2.set_xlim([np.min(LON), np.max(LON)])
ax2.set_ylim([np.min(LAT), np.max(LAT)])
plt.colorbar(im2)
plt.title('GLBb0.08 topo sub to IAS '+ftopo_sub)

dH = HH-HHglb
clrmp2 = plt.cm.Spectral
clrmp2.set_over(color=[0.3,0.3,0.3])
plt.figure(2,figsize=(8,8))
plt.clf()
ax3 = plt.subplot(2,1,1)
im3 = plt.pcolormesh(LON,LAT,dH,shading='flat',cmap=clrmp2,vmin=-1000., vmax=1000.)
plt.colorbar(im3)
plt.title('HH IAS - HH GLB')

heias = HH[:,-2]
heias[heias>100.]=np.nan
heglb = HHglb[:,-2]
heglb[heglb>100.]=np.nan

ax4 = plt.subplot(4,1,3)
ax4.plot(heias,label='Hias')
ax4.plot(heglb,label='Hglb')
ax4.legend()
plt.title('East OB')

heias = HH[-2,:]
heias[heias>100.]=np.nan
heglb = HHglb[-2,:]
heglb[heglb>100.]=np.nan
ax5 = plt.subplot(4,1,4)
ax5.plot(heias,label='Hias')
ax5.plot(heglb,label='Hglb')
ax5.legend()
plt.title('North OB')

btx = 'check_subGLBtopo_ias.py'
bottom_text(btx)







