"""
	Find indces for subsetting Global GOFS3.1 --> IAS domain
	1st step for creating nest flies
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

pthglb = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/misc/'
ftopo_glb = 'depth_GLBb0.08_09m11'
fgrid_glb = 'regional.grid'

import mod_read_hycom
from mod_read_hycom import read_grid_topo
lon_ias, lat_ias, hh_ias = read_grid_topo(pthias,ftopo_ias,fgrid_ias)
lon_glb, lat_glb, hh_glb = read_grid_topo(pthglb,ftopo_glb,fgrid_glb)

lon_min = np.min(lon_ias)
lon_max = np.max(lon_ias)
lat_min = np.min(lat_ias)
lat_max = np.max(lat_ias)

import mod_misc1
importlib.reload(mod_misc1)
from mod_misc1 import dist_sphcrd

# SW corner:
D = dist_sphcrd(lat_min,lon_min,lat_glb,lon_glb)
ng = D.shape[1]
mg = D.shape[0]
I = np.argmin(D)
ii = np.unravel_index(I,(mg,ng))
jg1 = ii[0]-1
ig1 = ii[1]-1

# NE corner:
D = dist_sphcrd(lat_max,lon_max,lat_glb,lon_glb)
I = np.argmin(D)
ii = np.unravel_index(I,(mg,ng))
jg2 = ii[0]+1
ig2 = ii[1]+1

print('IAS min lat={0}, min lon={1}'.format(lat_min,lon_min))
print(' Located outside closest Global pnts min lat={0}, min lon={1}'.\
			format(lat_glb[jg1,ig1],lon_glb[jg1,ig1]))
print('IAS max lat={0}, max lon={1}'.format(lat_max,lon_max))
print(' Located outside closest Global pnts max lat={0}, max lon={1}'.\
			format(lat_glb[jg2,ig2],lon_glb[jg2,ig2]))


#
# Subsample Glb topo and grid:
HH = hh_glb[jg1:jg2,ig1:ig2]
LON = lon_glb[jg1:jg2,ig1:ig2]
LAT = lat_glb[jg1:jg2,ig1:ig2]

# create figure and axes instances
Hp=HH.copy()

Hp[Hp>1.e20]=np.nan
Hp = -1.*Hp

f_plt = False
if f_plt:
	plt.ion()  # enables interactive mode
	fig1 = plt.figure(1,figsize=(8,8))
	plt.clf()
	ax1 = plt.subplot(1,1,1)
	im1 = plt.pcolormesh(LON,LAT,Hp,shading='flat',cmap='winter')
	ax1.axis('scaled')
	ax1.set_xlim([np.min(LON), np.max(LON)])
	ax1.set_ylim([np.min(LAT), np.max(LAT)])

	plt.colorbar(im1)
#
# Write indices
findx = pthglb+'indx_Glb2IAS.dat'
#findx_out = open(findx,'w')
ls = ['{0} iSW column index'.format(ig1), \
			'{0} jSW row index'.format(jg1), \
			'{0} iNE column index'.format(ig2), \
			'{0} jNE row index'.format(jg2)]

with open(findx,'w') as findx_out:
	for sls in ls[0:]:
		line = '{0}\n'.format(sls)
		findx_out.write(line)


# 
# Write modified topo for nest
# import struct
ftopo_sub = pthglb+'depth_GLBb0.08subIAS_09m11.a'
fgrid_sub = pthglb+'regional.grid_GLBb0.08subIAS.a'
print('Writing topo GLB sub IAS -->',ftopo_sub)
ftopo_out = open(ftopo_sub,'wb')
HH = HH.flatten()
IJDM = HH.shape[0]
# big endian, little: A.astype('f').byteswap().tofile('A.bin')
HH.astype('>f').tofile(ftopo_out) 
npad = 4096-IJDM%4096
toto = np.ones(npad,dtype='float32')
toto.astype('>f').tofile(ftopo_out)
ftopo_out.close()


print('Writing grid GLB sub IAS -->',fgrid_sub)
fgrid_out = open(fgrid_sub,'wb')
LON = LON.flatten()

LON.astype('>f').tofile(fgrid_out)
toto.astype('>f').tofile(fgrid_out)
LAT.astype('>f').tofile(fgrid_out)
toto.astype('>f').tofile(fgrid_out)

fgrid_out.close()









