"""
	Check oceanic vol transports at the OBs in the GLBb0.08 GOFS3.1 expt 73.7 
	used to create nest files for IAS HYCOM 0.03 
	2 approaches: (1) interpolate GLBb bathymetry to 0.03 to create nest topo
	and then interpolate archm fields to GLB-IAS0.03 topo
	In this approach, I would need to create "smoothed" bathymetry for IAS0.03
	to merge interior IAS topo to interpoalted form GLBb along the relaxation zone
	(2) interpolate archm fields directly onto IAS topo
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
from mod_read_hycom import read_hycom

def calc_vflux(U,dh,dx):
	"""
	Calcualte vol flux across the 2D section
	"""
	DX,dmm = np.meshgrid(dx,np.arange(1,ll+1))

	UFlx = U*DX*dh  # m3/s
	Utot = np.nansum(UFlx,axis=0)*1.e-6  # depth-integrated volume transport, Sv

	return Utot


def get_obs(F,iw=0,ie=-1,js=0,jn=-1):
	"""
	Extract fields at the OBS: N,S,E
	"""
	Aeast  = np.squeeze(F[:,ie,:])
	Asouth = np.squeeze(F[js,:,:])
	Anorth = np.squeeze(F[jn,:,:])

# Depth - first
	Aeast=Aeast.transpose()
	Asouth=Asouth.transpose()
	Anorth=Anorth.transpose()	

	return Aeast, Asouth, Anorth


# Fields from nest files on IAS0.03 topo:
LON, LAT, HH = read_grid_topo(pthias,ftopo_ias,fgrid_ias)
#HH[HH<1.e10] = -1.*HH[HH<1.e10]
mm = HH.shape[0]
nn = HH.shape[1]

pthbin = '/Net/kronos/ddmitry/hycom/TSIS/nest_files_gofs3.1/'
#fina = pthbin+'archm.ias003topo_2019_001.a'
#finb = pthbin+'archm.ias003topo_2019_001.b'

# Nest files from Navy HPC koehr same but different naming:
fina = pthbin+'archm.GLBb008ias003topo_001.a'
finb = pthbin+'archm.GLBb008ias003topo_001.b'

F, IDM, JDM, ll = read_hycom(fina,finb,'u-vel.')
Ueast, Usouth, Unorth = get_obs(F)

F, IDM, JDM, ll = read_hycom(fina,finb,'v-vel.')
Veast, Vsouth, Vnorth = get_obs(F)

# Layer thicknesses:
rg = 9806.  # convert pressure to depth, m
huge = 1.e20
F, IDM, JDM, ll = read_hycom(fina,finb,'thknss')
F[F>huge] = np.nan
F = F/rg
dHeast, dHsouth, dHnorth = get_obs(F)

import mod_utils_map
from mod_utils_map import dx_dy_2Dgrid
DX, DY = dx_dy_2Dgrid(LON,LAT)

DXeast  = DX[:,-1]
DYeast  = DY[:,-1]
DXsouth = DX[0,:]
DYsouth = DY[0,:]
DXnorth = DX[-1,:]
DYnorth = DY[-1,:]

# Flux from nest file with original IAS topo
UFe_ias = calc_vflux(Ueast,dHeast,DYeast)
UFs_ias = calc_vflux(Usouth,dHsouth,DXsouth)
UFn_ias = calc_vflux(Unorth,dHnorth,DXnorth)
He_ias = HH[:,-1]
Hs_ias = HH[0,:]
Hn_ias = HH[-1,:]

#
# Fields from the nest files with GLB0.08 09m11 topo
# interpolated onto IAS 0.03
LON, LAT, HH = read_grid_topo(pthias,ftopo_sub,fgrid_ias)
#HH[HH<1.e10] = -1.*HH[HH<1.e10]

fina = pthbin+'archm.topo09m11ias003_2019_001.a'
finb = pthbin+'archm.topo09m11ias003_2019_001.b'

F, IDM, JDM, ll = read_hycom(fina,finb,'u-vel.')
Ueast, Usouth, Unorth = get_obs(F)

F, IDM, JDM, ll = read_hycom(fina,finb,'v-vel.')
Veast, Vsouth, Vnorth = get_obs(F)

# Layer thicknesses:
rg = 9806.  # convert pressure to depth, m
huge = 1.e20
F, IDM, JDM, ll = read_hycom(fina,finb,'thknss')
F[F>huge] = np.nan
F = F/rg
dHeast, dHsouth, dHnorth = get_obs(F)

# Flux from nest file with GLBb topo--> 0.03IAS topo
UFe_glb = calc_vflux(Ueast,dHeast,DYeast)
UFs_glb = calc_vflux(Usouth,dHsouth,DXsouth)
UFn_glb = calc_vflux(Unorth,dHnorth,DXnorth)
# Chop off Pacific Ocean in nest with GLB topo:
UFs_glb[:1000]=np.nan
He_glb = HH[:,-1]
Hs_glb = HH[0,:]
Hn_glb = HH[-1,:]


plt.ion()  # enables interactive mode
fig = plt.figure(1, figsize=(8,8))
fig.clf()
ax = plt.axes([0.1, 0.7,0.8,0.25])
ax.plot(UFn_ias,label='IAStopo')
ax.plot(UFn_glb,label='09m11')
Fnias = np.nansum(UFn_ias)
Fnglb = np.nansum(UFn_glb)
plt.title('VolTrt, Sv, North, IAS={0:5.1f}, GLB={1:5.1f}, err={2:.3f}'.\
					format(Fnias,Fnglb,abs(Fnias-Fnglb)/abs(Fnglb)))
ax.legend(loc='upper left')

ax2 = plt.axes([0.1, 0.4,0.8,0.25])
ax2.plot(UFe_ias,label='IAStopo')
ax2.plot(UFe_glb,label='09m11')
Feias = np.nansum(UFe_ias)
Feglb = np.nansum(UFe_glb)
plt.title('VolTrt, Sv, North, IAS={0:5.1f}, GLB={1:5.1f}, err={2:.3f}'.\
          format(Feias,Feglb,abs(Feias-Feglb)/abs(Feglb)))

ax3 = plt.axes([0.1, 0.08,0.8,0.25])
ax3.plot(UFs_ias,label='IAStopo')
ax3.plot(UFs_glb,label='09m11')
plt.xlim([1300,nn])
Fsias = np.nansum(UFs_ias)
Fsglb = np.nansum(UFs_glb)
plt.title('VolTrt, Sv, North, IAS={0:5.1f}, GLB={1:5.1f}, err={2:.3f}'.\
          format(Fsias,Fsglb,abs(Fsias-Fsglb)/abs(Fsglb)))

btx = 'check_OBtransport.py'
bottom_text(btx, pos=[0.05,0.02])


f_chck = False
if f_chck:
	fig = plt.figure(1, figsize=(8,8))
	fig.clf()
	ax = plt.axes([0.1,0.3,0.8,0.6])
	im = plt.pcolormesh(F, shading='flat')
	plt.clim([-0.25, 0.25])
	ax.axis('scaled')

	plt.colorbar(im)

