# Imports
import numpy as np
import sys

# Paths to the .a and .b files
a = 'C:/Users/salfr/Downloads/GGS_TEST/rtofs_glo.t00z.f12.archv.a'
b = 'C:/Users/salfr/Downloads/GGS_TEST/rtofs_glo.t00z.f12.archv.b'

# Function to read the .a and .b files
def read_hycom(fina,finb,fld,Rtrc=None,rLayer=None,finfo=True):
  """
  reads hycom binary archive files (model output), 
  returns specified field 'fld'
  and dimensions of the grid
  Rtrc - passive tracer # to read, if there are tracers 
  rLayer - layer number to read, otherwise all layers are read - more time/memory
          numbering is lr=1,...,nlayers in the model
          rLayer is a layer number not an index, i.e.
          rLayer = 1 is 1st layer corresponds to index 0 in python

%  2017: added options for nlayer, n tracers
% if several tracers, options are:
% read tracer N: 'r_tracer',1
% read all tracers by default
% any variable can be read in 1 layer Nl:
%       'r_layer',1
%
% If all tracers are read, recommended to specify 
% only 1 layer to read 
%
  """
  try: 
    fgb = open(finb,'r')
  except:
    print('Could not open '+finb)

  fgb.seek(0)
  nl0 = 0
  while nl0 < 100:
    nl0 += 1
    data = fgb.readline().split()
    if len(data) < 2:
      continue
    adim = data[1]
    ii = adim.find('idm')
    if ii>0:
      break

  if ii<0:
    fgb.close()
    sys.exit('No idm found: Reading ',finb)

  IDM = int(data[0])
  data = fgb.readline().split()
  JDM = int(data[0])
  IJDM = IDM*JDM
#  fgb.seek(0)

  npad =4096-IJDM%4096

  if finfo:
    print('Reading HYCOM :{0} '.format(finb))

  aa = fgb.readline().split()

# Find location of the requested field:
  cntr= 0
  nf = len(fld)
  FLOC=[]
  while cntr<1e6:
    aa = fgb.readline().split()
    if len(aa) == 0:  # end of file
      break
    cntr += 1
    aname = aa[0]
    ii = aname.find(fld)
    if ii >= 0 and len(aname)==nf:
      FLOC.append(cntr)

  fgb.close()

  nrec = len(FLOC)
  if nrec == 0:
    raise Exception('read_hycom: Field {0} not found in {1}'.format(fld,finb))

# N. of v. layers
  """
   If fld = tracer and # of tracers >1
   need to distinguish # of layers 
   vs # of tracers*layers
   if strmatch(fld,'tracer')
  """
  ll = len(FLOC)
  if finfo:
    print('Grid: IDM={0}, JDM={1}, KDM={2}'.format(IDM,JDM,ll))

  FLOC = np.array(FLOC)
  if ll == 1:
    nVlev = 1
    nTR = 1
  else:
    dI = np.diff(FLOC)
    dmm = np.where(dI>1)
    dindx = dmm[0]
    nTR = dindx[0]+1         # N. of tracers in 1 layer
    nVlev = ll/nTR        # # of v. layers

# breakpoint()
  if nTR != 1 and finfo:
    print('read_hycom: Found {0} variables {1}  per layer'.format(nTR,fld))

# Find layers to read, if specified
# and tracers, if specified
  lr1=-1
  lr2=-1
  if Rtrc is not None:
    if nTR < Rtrc:
      raise Exception('Number of saved tracers {0} < requested {1}'.\
                      format(nTR,Rtrc))
    dmm = np.copy(FLOC)
    FLOC = dmm[Rtrc-1::nTR]

    if lr1 < 0 or lr2 < 0 :
      lr2 = FLOC.shape[0]
      ll = lr2

  if rLayer is not None:
    lr1 = rLayer
    lr2 = lr1

# If a layer Number not specified - read all
  if lr1 < 0 or lr2 < 0:
    lr1 = 1
    lr2 = ll

#  print('Reading {0}, Layers: {1}-{2}'.format(fld,lr1,lr2))
  fga = open(fina,'rb')
  F = []
  ccL = -1
  huge = 0.0001*2.**99
  for ii in range(lr1,lr2+1):
    fga.seek(0)
    k0 = FLOC[ii-1]-1
    fga.seek(np.int64(k0) * (npad + IJDM) * 4,0)
    dmm = np.fromfile(fga, dtype='>f',count=IJDM) # read 1 layer
    amin = np.min(dmm[np.where(dmm<huge)])
    amax = np.max(dmm[np.where(dmm<huge)])
    if finfo:
      print('Reading {0} k={1} min={2:12.4f} max={3:12.4f}'.\
            format(fld,ii,amin,amax))
    dmm = dmm.reshape((JDM,IDM))
    ccL += 1
#   print('lr={0}, ccL={1}'.format(ii,ccL))
#   breakpoint()
    if ccL == 0:
      F = np.copy(dmm)
    else:
      F = np.dstack((F,dmm))

  if ll == 0:
    print('!!! read_hycom: {0} not found in {1} ERR'.format(fld,fina))
    print('!!! read hycom: check fields in ',finb)

  fga.close()

  return F, IDM, JDM, ll

####################################################################################################

import matplotlib.pyplot as plt
import numpy as np

def sub_check_sect(F, Fi, ZM, ZMi, Hs):
    """
    Plot 2 sections for checking
    """
    
    ll, msb = F.shape
    xx, dmm = np.meshgrid(range(1, msb+1), range(1, ll+1))
    
    # Plot for F
    plt.figure(11)
    plt.clf()
    plt.pcolormesh(xx, ZM, F, shading='flat')
    plt.colorbar()
    plt.clim(np.min(F), np.max(F))

    ll, msb = Fi.shape
    xxi, zzi = np.meshgrid(range(1, msb+1), ZMi)
    
    # Plot for Fi
    plt.figure(12)
    plt.clf()
    plt.pcolormesh(xxi, zzi, Fi, shading='flat')
    plt.hold(True)  # This is not usually required in newer versions of matplotlib
    plt.plot(xxi, Hs, 'r-')
    plt.colorbar()
    plt.clim(np.min(F), np.max(F))
    plt.title('Interpolated into z-layers')

    plt.show()

####################################################################################################

import numpy as np

def sub_fill_3z(F):
    """
    3D array
    Fill all nans below the 1st upper layer
    with values for vertical interpolation
    1 layer must not have nans
    """
    print("Filling nans 3z")

    # Check if the first layer has NaN values
    a = F[0, :, :]
    if np.isnan(a).any():
        raise ValueError("sub_fill_3z: 1st layer cannot have nans")

    # Get the shape of F
    kk, mm, nn = F.shape

    # Fill NaNs
    for ik in range(1, kk):
        In = np.isnan(F[ik, :, :])
        F[ik][In] = F[ik-1][In]

    return F

####################################################################################################

import numpy as np

def sub_fill_bottom_nans(Ssm):
    """
    2D array
    For plotting and interpolation
    Fill bottom nans with closest values
    """
    # Check if the entire array has no NaN values
    if not np.isnan(Ssm).any():
        return Ssm

    print("Filling bottom nans ...")

    mm, nn = Ssm.shape

    # Check and fill NaNs on the left border
    if np.isnan(Ssm[0, 0]):
        i1 = np.min(np.where(~np.isnan(Ssm[0, :])))
        for kk in range(mm):
            Ssm[kk, :i1] = Ssm[kk, i1]

    # Check and fill NaNs on the right border
    if np.isnan(Ssm[-1, 0]):
        i1 = np.max(np.where(~np.isnan(Ssm[0, :])))
        for kk in range(mm):
            Ssm[kk, i1+1:] = Ssm[kk, i1]

    # Compute the mean ignoring NaNs
    s0 = np.nanmean(Ssm)

    for ii in range(nn):
        i1_values = np.where(np.isnan(Ssm[:, ii]))[0]
        if i1_values.size > 0:
            i1 = i1_values[0]
            if i1 > 0:
                Ssm[i1:, ii] = Ssm[i1-1, ii]
            elif i1 == 0:  # island
                Ssm[:, ii] = Ssm[:, ii-1]

    return Ssm

####################################################################################################

import numpy as np

def sub_fill_land(F, IFL=None):
    """
    Fill land areas in the horizontal plane using closest ocean values.

    Parameters:
        F (ndarray): Input 2D array.
        IFL (dict, optional): Indices data structure. If not provided, it will be computed.

    Returns:
        tuple: Modified F array and IFL dictionary.
    """
    ma, na = F.shape

    if IFL is None:
        IFL = {}
        In = np.argwhere(np.isnan(F))
        JJ, II = np.argwhere(~np.isnan(F)).T

        print("Deriving Fill Land indices ...")

        IFL["Land_Indx"] = []
        IFL["Land_IJ"] = []
        IFL["Ocean_Indx"] = []
        IFL["Ocean_IJ"] = []

        for ik, (j0, i0) in enumerate(In):
            if ik % 10000 == 0:
                print(f"  Done {ik / len(In) * 100:.2f}% ...")

            d = np.sqrt((JJ - j0)**2 + (II - i0)**2)
            kk = np.argmin(d)
            jm, im = JJ[kk], II[kk]

            IFL["Land_Indx"].append(ik)
            IFL["Land_IJ"].append((i0, j0))
            IFL["Ocean_Indx"].append(np.ravel_multi_index((jm, im), (ma, na)))
            IFL["Ocean_IJ"].append((im, jm))

        # Convert lists to arrays for better indexing performance
        IFL["Land_Indx"] = np.array(IFL["Land_Indx"])
        IFL["Ocean_Indx"] = np.array(IFL["Ocean_Indx"])

    # Fill land areas with ocean values
    F[IFL["Land_Indx"]] = F[IFL["Ocean_Indx"]]

    return F, IFL

####################################################################################################

import numpy as np
from scipy.interpolate import pchip_interpolate

def sub_interp2z_2D(Hs, F, ZMf, ZMh):
    """
    Interpolate normal to a vertical section (2D array) from HYCOM, hybrid grid onto z-fixed depths.
    
    Parameters:
        Hs (ndarray): Bottom depth.
        F (ndarray): Input 2D field.
        ZMf (ndarray): 1D array of mid-layer depths for fixed layers.
        ZMh (ndarray): HYCOM 2D array of rho-point depths.
        
    Returns:
        ndarray: Interpolated field Fi.
    """
    hg = 1e20
    nlr, nn = F.shape
    mZ = len(ZMf)
    np = len(Hs)

    # Add extra bottom layer for interpolation below the min depth on fixed Z-array
    ZMh = np.vstack((ZMh, -12000 * np.ones((1, nn))))
    F = np.vstack((F, F[-1, :]))

    # Add extra surface layer for interpolation
    ZMh = np.vstack((np.zeros((1, nn)), ZMh))
    F = np.vstack((F[0, :], F))

    # If there are NaNs in F, use sub_fill_bottom_nans function (to be defined separately)
    # NOTE: You must ensure the function `sub_fill_bottom_nans` is defined in your Python code.
    if np.isnan(F).any():
        F = sub_fill_bottom_nans(F)

    nlr, nn = F.shape
    Xh, dmm = np.meshgrid(np.arange(1, nn + 1), np.arange(1, nlr + 1))

    # Interpolation
    Fi = np.zeros((mZ, np))
    for kk in range(np):
        zm = ZMh[:, kk]
        
        # For interpolation, depths cannot be the same
        dz = np.abs(np.diff(zm))
        Iz = np.where(dz == 0)[0]
        for iz0 in Iz:
            zm[iz0 + 1] = zm[iz0] - 0.001
        
        f0 = F[:, kk]
        fi = pchip_interpolate(zm, f0, ZMf)
        Fi[:, kk] = fi

    return Fi

####################################################################################################

import numpy as np
from scipy.interpolate import interp2

def sub_interp3D_nas(T, IFL, XX, YY, HHs, ZMnas, ZM, nlon, nlat, ILAND, pfld):
    """
    Interpolate 3D field (T) onto NAS grid.
    """
    huge = 1e30

    msb, nsb = XX.shape
    knas = len(ZMnas)
    lnas = len(ZMnas)
    mnas, nnas = nlon.shape

    # fill NaNs in the surface layer
    t1 = T[0, :, :]
    t1, IFL = sub_fill_land(t1, IFL)  # You must ensure `sub_fill_land` is defined in your Python code
    T[0, :, :] = t1

    # Fill NaNs in the subsurface layers (bottom)
    T = sub_fill_3z(T)  # You must ensure `sub_fill_3z` is defined in your Python code

    # Interpolate into z levels as a 2D vertical section
    Tzi = np.zeros((len(ZMnas), T.shape[1], nsb))

    for ii in range(nsb):
        if ii % 200 == 0:
            print(f"{pfld} interp to z, done {(ii / nsb) * 100:.2f}% ...")

        y1 = YY[:, ii]
        ZMh = ZM[:, :, ii]
        F = T[:, :, ii]
        Hs = HHs[:, ii]
        Hs[Hs >= 0] = -1

        Fi = sub_interp2z_2D(Hs, F, ZMnas, ZMh)  # You must ensure `sub_interp2z_2D` is defined

        # If you need to check section, uncomment below
        # sub_check_sect(F, Fi, ZMh, ZMnas, Hs)  # You must ensure `sub_check_sect` is defined

        Tzi[:, :, ii] = Fi

    # Interpolate onto NAS grid
    print(f"{pfld} Interpolating onto XY grid")
    a1, a2, a3 = Tzi.shape
    Tzx = np.zeros((a1, mnas, nnas))

    for kk in range(a1):
        F = Tzi[kk, :, :]
        Fi = interp2(XX, YY, F, nlon, nlat)
        Tzx[kk, :, :] = Fi

    # Put land back
    for kk in range(knas):
        il = ILAND[kk]['I']
        Tzx[kk, il] = huge

    return Tzx

####################################################################################################

import netCDF4 as nc
from datetime import datetime

def sub_write_netcdf(flcdf, nlon, nlat, Znas, sshi, Tzx, Szx, Uzx, Vzx, dnmb, huge):
    """
    Create and write data to a new netCDF file.
    """
    
    # Calculate days since 1900-12-31
    d0 = datetime(1900, 12, 31)
    delta = dnmb - d0
    nday = delta.days

    mnas, nnas = nlon.shape
    knas = len(Znas)

    print(f'Writing {flcdf}')

    # Open netCDF file for writing
    ncid = nc.Dataset(flcdf, 'w', format='NETCDF4')

    # Define dimensions
    dimidt = ncid.createDimension('time', 1)
    dimidy = ncid.createDimension('north', mnas)
    dimidx = ncid.createDimension('east', nnas)
    dimidz = ncid.createDimension('depth', knas)

    # Define variables
    dayID = ncid.createVariable('time', 'f8', ('time',))
    lonID = ncid.createVariable('longitude', 'f8', ('east', 'north'))
    latID = ncid.createVariable('latitude', 'f8', ('east', 'north'))
    dpthID = ncid.createVariable('depth', 'f8', ('depth',))
    sshID = ncid.createVariable('srfhgt', 'f8', ('east', 'north'))
    tempID = ncid.createVariable('temp', 'f8', ('east', 'north', 'depth'))
    salnID = ncid.createVariable('saln', 'f8', ('east', 'north', 'depth'))
    uvelID = ncid.createVariable('u-vel', 'f8', ('east', 'north', 'depth'))
    vvelID = ncid.createVariable('v-vel', 'f8', ('east', 'north', 'depth'))

    # Set Fill Value
    sshID.setncattr('_FillValue', huge)
    tempID.setncattr('_FillValue', huge)
    salnID.setncattr('_FillValue', huge)
    uvelID.setncattr('_FillValue', huge)
    vvelID.setncattr('_FillValue', huge)

    # Global Attributes
    ncid.info = 'HYCOM-TSIS GoM 0.03, COAPS FSU'
    ncid.creation_date = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # Write Time Data
    dayID[:] = nday
    dayID.standard_name = 'Days'
    dayID.units = 'days since 1900-12-31 00:00:00'

    # Write Grid Data
    lonID[:] = nlon.T
    lonID.standard_name = 'Longitude'
    lonID.units = 'degrees_east'
    lonID.axis = 'X'

    latID[:] = nlat.T
    latID.standard_name = 'Latitude'
    latID.units = 'degrees_north'
    latID.axis = 'Y'

    dpthID[:] = Znas
    dpthID.standard_name = 'Depths of fixed Z levels'
    dpthID.units = 'm'
    dpthID.axis = 'Z'

    # Write SSH
    sshID[:] = sshi.T
    sshID.standard_name = 'sea surface height'
    sshID.units = 'm'

    # Write 3D Data
    tempID[:] = np.transpose(Tzx, (2, 1, 0))
    tempID.standard_name = 'potential temperature'
    tempID.units = 'degrees Celsius'

    salnID[:] = np.transpose(Szx, (2, 1, 0))
    salnID.standard_name = 'salinity'
    salnID.units = 'salinity units'

    uvelID[:] = np.transpose(Uzx, (2, 1, 0))
    uvelID.standard_name = 'u-component of ocean velocity vector'
    uvelID.units = 'm/s'
    uvelID.positive_dir = 'eastward'

    vvelID[:] = np.transpose(Vzx, (2, 1, 0))
    vvelID.standard_name = 'v-component of ocean velocity vector'
    vvelID.units = 'm/s'
    vvelID.positive_dir = 'northward'

    ncid.close()

    print('File created')

####################################################################################################

import mat2py as mp
from mat2py.core import *

def main():
    addpath("/usr/people/ddmitry/codes/MyMatlab/")
    addpath("/usr/people/ddmitry/codes/MyMatlab/hycom_utils")
    addpath("/usr/people/ddmitry/codes/MyMatlab/colormaps")
    startup
    close("all")
    clear
    f_write = 1
    expt = "noPIES"
    ys = 2011
    ye = 2011
    rg = 9806
    huge = 1e25
    fprintf(" %s %i-%i\n", expt, ys, ye)
    pthfig = "/Net/mars/ddmitry/hycom/hycom_TSIS/fig_ssh/"
    pthtopo = "/home/ddmitry/codes/HYCOM_TSIS/"
    pthgrid = "/Net/mars/ddmitry/hycom/hycom_TSIS/data_project/"
    btx = "interp_hycom2nas_grid.m"
    ftopo = sprintf("%sias_gridinfo.nc", pthtopo)
    HH = (-1) * nc_varget(ftopo, "mdepth")
    LAT = nc_varget(ftopo, "mplat")
    LON = nc_varget(ftopo, "mplon")
    mm, nn = size(HH)
    m = copy(mm)
    n = copy(nn)
    HH[I[isnan(HH)]] = 100
    DX, DY = sub_dx_dy(LON, LAT)
    XM, YM = meshgrid(M[M[1:n]], M[M[1:m]])
    fgrid = sprintf("%scommon_gom_domain.mat", pthgrid)
    NAS = load(fgrid)
    nlat = NAS.lat
    nlon = NAS.lon
    ynas1 = min(min(nlat))
    ynas2 = max(max(nlat))
    xnas1 = min(min(nlon))
    xnas2 = max(max(nlon))
    ZMnas = -(NAS.z)
    J, _I = find(LON < xnas1)
    its1 = max(_I)
    J, _I = find(LON > xnas2)
    its2 = min(_I)
    J, _I = find(LAT < ynas1)
    jts1 = max(J)
    J, _I = find(LAT > ynas2)
    jts2 = min(J)
    xt1 = LON(jts1, its1)
    yt1 = LAT(jts1, its1)
    xt2 = LON(jts2, its2)
    yt2 = LAT(jts2, its2)
    XX = LON(M[jts1:jts2], M[its1:its2])
    YY = LAT(M[jts1:jts2], M[its1:its2])
    HHs = HH(M[jts1:jts2], M[its1:its2])
    msb, nsb = size(HHs)
    knas = length(ZMnas)
    lnas = length(ZMnas)
    mnas, nnas = size(nlon)
    p_map = 0
    if p_map == 1:
        figure(10)
        clf
        contour(LON, LAT, HH, M[[0, 0]], "k")
        hold("on")
        plot(M[[xnas1, xnas1]], M[[ynas1, ynas2]], "r")
        plot(M[[xnas2, xnas2]], M[[ynas1, ynas2]], "r")
        plot(M[[xnas1, xnas2]], M[[ynas1, ynas1]], "r")
        plot(M[[xnas1, xnas2]], M[[ynas2, ynas2]], "r")
        plot(M[[xt1, xt1]], M[[yt1, yt2]], "g")
        plot(M[[xt2, xt2]], M[[yt1, yt2]], "g")
        plot(M[[xt1, xt2]], M[[yt1, yt1]], "g")
        plot(M[[xt1, xt2]], M[[yt2, yt2]], "g")
    for year in M[ys:ye]:
        pthd = sprintf(
            "/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/%i_%s/", year, expt
        )
        pthout = sprintf(
            "/nexsan/people/ddmitry/GOM_TSIS0.03/hindcast/%s/%4.4i/", expt, year
        )
        if _not(exist(pthout, "dir")):
            system(sprintf("mkdir -pv %s", pthout))
        id1 = 1
        id2 = 365
        if (mod(year, 4) == 0) & (id2 >= 365):
            id2 = 366
        for iday in M[id1:id2]:
            sday = sprintf("%3.3i", iday)
            hr = 0
            fina = sprintf("%sarchv.%4.4i_%3.3i_%2.2i.a", pthd, year, iday, hr)
            finb = sprintf("%sarchv.%4.4i_%3.3i_%2.2i.b", pthd, year, iday, hr)
            fin = copy(fina)
            ie = exist(fin, "file")
            if _not(ie):
                fprintf("Missing: %s\n", fin)
                continue
            tic
            dJ1 = datenum(year, 1, 1)
            nday = (dJ1 + iday) - 1
            dnmb = nday + (hr / 24)
            DV = datevec(dnmb)
            fprintf("Reading %s\n", fina)
            fprintf("\n iday=%i, %4.4i/%2.2i/%2.2i:%2.2ih\n", iday, DV(M[1:4]))
            date_str = sprintf("%4.4i/%2.2i/%2.2i : %2.2ih", DV(M[1:4]))
            if _not(exist("Lmask", "var")):
                LL = HHs * 0
                LL[I[HHs < 0]] = 1
                Ioc = find(LL == 1)
                Ild = find(LL == 0)
                Lmask = interp2(XX, YY, LL, nlon, nlat)
                Lmask = round(Lmask)
                HHs[I[Ild]] = 0.0
                Hnas = interp2(XX, YY, HHs, nlon, nlat)
                HHs[I[Ild]] = 100
                Hnas[I[Lmask == 0]] = 100
                Ib = find((Hnas >= 0) & (Lmask == 1))
                if _not(isempty(Ib)):
                    Hnas[I[Ib]] = -1
                for kk in M[1:knas]:
                    _I = find(Hnas >= ZMnas(kk))
                    ILAND(kk)._I = copy(_I)

            F, n, m, l = read_hycom(fina, finb, "thknss")
            F = mrdivide(F, rg)
            F[I[F > 1e10]] = 0
            dH = squeeze(F)
            ih = find(dH < 1e-1)
            dHs = dH[I[:, jts1:jts2, its1:its2]]
            ZM1, ZZ1 = sub_zz_zm(fina, finb, HH, "thknss", dH)
            ZM1[I[isnan(ZM1)]] = -100
            ZZ = ZZ1[I[:, jts1:jts2, its1:its2]]
            ZM = ZM1[I[:, jts1:jts2, its1:its2]]
            fld = "srfhgt"
            F, nn, mm, ll = read_hycom(fina, finb, fld)
            F[I[F > huge]] = copy(nan)
            ssh = squeeze(F) / (1e-3 * rg)
            ssh = ssh(M[jts1:jts2], M[its1:its2])
            if _not(exist("IFL", "var")):
                IFL = M[[]]
            fprintf("Interpolating ssh \n")
            ssh, IFL = sub_fill_land(ssh, IFL)
            sshi = interp2(XX, YY, ssh, nlon, nlat)
            il = ILAND(1)._I
            sshi[I[il]] = copy(huge)
            fld = "temp"
            F, nn, mm, ll = read_hycom(fina, finb, fld)
            F[I[F > huge]] = copy(nan)
            F = squeeze(F)
            T = F[I[:, jts1:jts2, its1:its2]]
            Tzx = sub_interp3D_nas(
                T, IFL, XX, YY, HHs, ZMnas, ZM, nlon, nlat, ILAND, fld
            )
            f_sct = 0
            if f_sct == 1:
                ii = 400
                ZMh = squeeze(ZM[I[:, :, ii]])
                F1 = squeeze(T[I[:, :, ii]])
                x1 = XX(1, ii)
                dd = abs(nlon[I[1, :]] - x1)
                ii2 = find(dd == min(dd), 1)
                F2 = squeeze(Tzx[I[:, :, ii2]])
                F2[I[F2 > 1e20]] = copy(nan)
                Hs = Hnas[I[:, ii2]]
                sub_check_sect(F1, F2, ZMh, ZMnas, Hs)
                bottom_text(btx, "pwd", 1)
            fld = "salin"
            F, nn, mm, ll = read_hycom(fina, finb, fld)
            F[I[F > huge]] = copy(nan)
            F = squeeze(F)
            S = F[I[:, jts1:jts2, its1:its2]]
            Szx = sub_interp3D_nas(
                S, IFL, XX, YY, HHs, ZMnas, ZM, nlon, nlat, ILAND, fld
            )
            f_sct = 0
            if f_sct == 1:
                ii = 400
                ZMh = squeeze(ZM[I[:, :, ii]])
                F1 = squeeze(S[I[:, :, ii]])
                x1 = XX(1, ii)
                dd = abs(nlon[I[1, :]] - x1)
                ii2 = find(dd == min(dd), 1)
                F2 = squeeze(Szx[I[:, :, ii2]])
                F2[I[F2 > 1e20]] = copy(nan)
                Hs = Hnas[I[:, ii2]]
                sub_check_sect(F1, F2, ZMh, ZMnas, Hs)
                bottom_text(btx, "pwd", 1)
            F, n, m = read_hycom(fina, finb, "u_btrop")
            F[I[F > 1e6]] = copy(nan)
            F = squeeze(F)
            Ub = F(M[jts1:jts2], M[its1:its2])
            fld = "u-vel"
            F, nn, mm, ll = read_hycom(fina, finb, fld)
            F[I[F > huge]] = copy(nan)
            F = squeeze(F)
            U = F[I[:, jts1:jts2, its1:its2]]
            for kk in M[1:ll]:
                U[I[kk, :, :]] = squeeze(U[I[kk, :, :]]) + Ub

            fprintf(" Collocating U\n")
            U[I[isnan(U)]] = 0
            Uclc = copy(U)
            for ii in M[2 : (nsb - 1)]:
                dh1 = squeeze(dHs[I[:, :, ii - 1]])
                dh2 = squeeze(dHs[I[:, :, ii]])
                dh12 = 0.5 * (dh1 + dh2)
                dh2 = squeeze(dHs[I[:, :, ii]])
                dh3 = squeeze(dHs[I[:, :, ii + 1]])
                dh23 = 0.5 * (dh1 + dh2)
                u12 = squeeze(U[I[:, :, ii]])
                u23 = squeeze(U[I[:, :, ii + 1]])
                u0 = ((u12 * dh12) + (u23 * dh23)) / (dh12 + dh23)
                a = dh12 + dh23
                Im = find(a < 1e-3)
                u0[I[Im]] = 0
                Uclc[I[:, :, ii]] = copy(u0)

            Uzx = sub_interp3D_nas(
                Uclc, IFL, XX, YY, HHs, ZMnas, ZM, nlon, nlat, ILAND, fld
            )
            f_sct = 0
            if f_sct == 1:
                ii = 400
                ZMh = squeeze(ZM[I[:, :, ii]])
                F1 = squeeze(U[I[:, :, ii]])
                x1 = XX(1, ii)
                dd = abs(nlon[I[1, :]] - x1)
                ii2 = find(dd == min(dd), 1)
                F2 = squeeze(Uzx[I[:, :, ii2]])
                F2[I[F2 > 1e20]] = copy(nan)
                Hs = Hnas[I[:, ii2]]
                sub_check_sect(F1, F2, ZMh, ZMnas, Hs)
                bottom_text(btx, "pwd", 1)
            F, n, m = read_hycom(fina, finb, "v_btrop")
            F[I[F > 1e6]] = copy(nan)
            F = squeeze(F)
            Vb = F(M[jts1:jts2], M[its1:its2])
            fld = "v-vel"
            F, nn, mm, ll = read_hycom(fina, finb, fld)
            F[I[F > huge]] = copy(nan)
            F = squeeze(F)
            V = F[I[:, jts1:jts2, its1:its2]]
            for kk in M[1:ll]:
                V[I[kk, :, :]] = squeeze(V[I[kk, :, :]]) + Vb

            fprintf(" Collocating V\n")
            V[I[isnan(V)]] = 0
            Vclc = copy(V)
            for jj in M[2 : (msb - 1)]:
                dh1 = squeeze(dHs[I[:, jj - 1, :]])
                dh2 = squeeze(dHs[I[:, jj, :]])
                dh12 = 0.5 * (dh1 + dh2)
                dh2 = squeeze(dHs[I[:, jj, :]])
                dh3 = squeeze(dHs[I[:, jj + 1, :]])
                dh23 = 0.5 * (dh1 + dh2)
                v12 = squeeze(V[I[:, jj, :]])
                v23 = squeeze(V[I[:, jj + 1, :]])
                v0 = ((v12 * dh12) + (v23 * dh23)) / (dh12 + dh23)
                a = dh12 + dh23
                Im = find(a < 1e-3)
                v0[I[Im]] = 0
                Vclc[I[:, jj, :]] = copy(v0)

            Vzx = sub_interp3D_nas(
                Vclc, IFL, XX, YY, HHs, ZMnas, ZM, nlon, nlat, ILAND, fld
            )
            f_sct = 0
            if f_sct == 1:
                ii = 400
                ZMh = squeeze(ZM[I[:, :, ii]])
                F1 = squeeze(V[I[:, :, ii]])
                x1 = XX(1, ii)
                dd = abs(nlon[I[1, :]] - x1)
                ii2 = find(dd == min(dd), 1)
                F2 = squeeze(Vzx[I[:, :, ii2]])
                F2[I[F2 > 1e20]] = copy(nan)
                Hs = Hnas[I[:, ii2]]
                sub_check_sect(F1, F2, ZMh, ZMnas, Hs)
                bottom_text(btx, "pwd", 1)
            if f_write == 1:
                flcdf = sprintf(
                    "%shycom0.03_%s_%4.4i_%3.3i.nc", pthout, expt, year, iday
                )
                sub_write_netcdf(
                    flcdf, nlon, nlat, ZMnas, sshi, Tzx, Szx, Uzx, Vzx, dnmb, huge
                )
            fprintf("1 day processed %6.2f min\n\n", toc / 60)

####################################################################################################



####################################################################################################