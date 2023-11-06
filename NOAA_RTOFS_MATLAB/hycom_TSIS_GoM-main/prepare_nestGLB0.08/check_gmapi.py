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

pthias = '/Net/kronos/ddmitry/hycom/TSIS/topo_IAS0.03/'
ftopo_ias = 'regional.depth'
fgrid_ias = 'regional.grid'
fgmapi = 'regional.gmapi_GLB0.08'

def read_gmapi(pthtopo,fgmapi,IDM,JDM):
	flgmapia = pthtopo+fgmapi+'.a'
	flgmapib = pthtopo+fgmapi+'.b'
	IJDM = IDM*JDM
  npad =4096-IJDM%4096
#
# Direct access:
	fid = open(flgmapia,'rb')
	fid.seek(0)
	x_out = np.fromfile(fid,dtype='>f',count=IJDM)
	fid.seek(4*(npad+IJDM),0)
	y_out = np.fromfile(fid,dtype='>f',count=IJDM)
	fid.close()

	x_out = x_out.reshape((JDM,IDM))
	y_out = y_out.reshape((JDM,IDM))

	print('gmapi: x_out min/max {0}, {1}'.format(np.min(x_out),np.max(x_out)))
	print('gmapi: y_out min/max {0}, {1}'.format(np.min(y_out),np.max(y_out)))

	return x_out, y_out

import mod_read_hycom
from mod_read_hycom import read_grid_topo
LON, LAT, HH = read_grid_topo(pthias,ftopo_ias,fgrid_ias)
mm = HH.shape[0]
nn = HH.shape[1]

x_out, y_out = read_gmapi(pthias,fgmapi,nn,mm)

plt.ion()  # enables interactive mode
fig1 = plt.figure(1,figsize=(8,8))
plt.clf()
ax1 = plt.subplot(2,1,1)
im1 = plt.pcolormesh(LON,LAT,x_out,shading='flat',cmap='winter')
ax1.axis('scaled')
ax1.set_xlim([np.min(LON), np.max(LON)])
ax1.set_ylim([np.min(LAT), np.max(LAT)])
plt.colorbar(im1)
plt.title('x_out '+fgmapi)

ax2 = plt.subplot(2,1,2)
im2 = plt.pcolormesh(LON,LAT,y_out,shading='flat',cmap='winter')
ax2.axis('scaled')
ax2.set_xlim([np.min(LON), np.max(LON)])
ax2.set_ylim([np.min(LAT), np.max(LAT)])
plt.colorbar(im2)
plt.title('y_out '+fgmapi)




