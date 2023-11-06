"""
	Prepare nest bathymetry file from the topo-file 
	IAS 0.03 by removing "walls" at the OBs
"""
from IPython import get_ipython
get_ipython().magic('reset -sf')

import os
from netCDF4 import Dataset as ncFile
import numpy as np
from copy import copy
import importlib
import matplotlib.pyplot as plt


drnm = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/restart/'
rname_out = 'cice.restart_117f-eap.nc'

pthtopo = '/Net/kronos/ddmitry/hycom/TSIS/topo_IAS0.03/'
pthnc = '/Net/kronos/ddmitry/hycom/TSIS/' 
ftoponc = 'ias_gridinfo.nc'
ftopo = 'regional.depth'
fgrid = 'regional.grid'


def get_grid(pthgrid,ftopo):
#  pthgrid = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/'
#  flgrd = 'depth_ARCc0.04_17DD.nc'
	nc = ncFile(pthgrid+ftopo)
	lonc = nc.variables['mplon'][:]
	latc = nc.variables['mplat'][:]
	lonc[lonc<0]=lonc[lonc<0]+360
	HH = nc.variables['Bathymetry'][:]
	HH = -1.*HH

	return latc,lonc,HH

def get_grid_ab(pthtopo,ftopo,fgrid):
# read HYCOM grid and topo files *.[ab]
	fltopoa = pthtopo+ftopo+'.a'
	fltopob = pthtopo+ftopo+'.b'
	flgrida = pthtopo+fgrid+'.a'
	flgridb = pthtopo+fgrid+'.b'
	
# Read dim, lon, lat
	fgb = open(flgridb,'r')
	fgb.seek(0)
	data = fgb.readline().split()
	IDM = int(data[0])
	data = fgb.readline().split()
	JDM = int(data[0])
	IJDM = IDM*JDM
	fgb.close()

	npad =4096-IJDM%4096

	print('Reading HYCOM grid/topo:{0} {1} ',fgrid,ftopo)
	print('Grid: IDM={0}, JDM={1}'.format(IDM,JDM))

# Read direct access binary file
	fga = open(flgrida,'rb')
	fga.seek(0)
	plon = np.fromfile(fga,dtype='>f',count=IJDM)
	I = np.where(plon > 180.)
	plon[I] = plon[I]-360.
	plon = plon.reshape((JDM,IDM))

	fga.seek(4*(npad+IJDM),0)
	plat = np.fromfile(fga, dtype='>f',count=IJDM)
	plat = plat.reshape((JDM,IDM))

	fga.close()

# Read bottom topography:
	fbt = open(fltopoa,'rb')
	fbt.seek(0)
	HH = np.fromfile(fbt, dtype='>f', count=IJDM)
	HH = HH.reshape((JDM,IDM))
	fbt.close()

	return plon,plat,HH

LON, LAT, HH = get_grid_ab(pthtopo,ftopo,fgrid)


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
# Add values at the closed OB
HH[:,-1]=HH[:,-2]
HH[-1,:-1]=HH[-2,:-1]
HH[-1,-1]=HH[-2,-1]

# Write modified topo for nest
#import struct
fltopo_nesta = pthtopo+'regional.depth-nest.a'
print('Writing nest topo -->',fltopo_nesta)
ftopo_out = open(fltopo_nesta,'wb')
HH = HH.flatten()
IJDM = HH.shape[0]
# big endian, little: A.astype('f').byteswap().tofile('A.bin')
HH.astype('>f').tofile(ftopo_out) 
npad = 4096-IJDM%4096
toto = np.ones(npad,dtype='float32')
toto.astype('>f').tofile(ftopo_out)
ftopo_out.close()







