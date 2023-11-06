import numpy as np
import sys

from sup_functions import dist_sphcrd

a = 'C:/Users/salfr/Downloads/GGS_TEST/rtofs_glo.t00z.f12.archv.a'
b = 'C:/Users/salfr/Downloads/GGS_TEST/rtofs_glo.t00z.f12.archv.b'

def read_ab(fina,finb,fld,Rtrc=None,rLayer=None,finfo=True):
    '''
    Reads hycom binary archive files (RTOFS/HYCOM model output)

    Args:
    - fina (str): The path to the .a file.
    - finb (str): The path to the .b file.
    - fld (str): The field to be read.
        - 'u-vel.': u-velocity
        - 'v-vel.': v-velocity
        - 'thknss': layer thickness
        - 'temp': temperature
        - 'salin': salinity
    - Rtrc (int): The tracer number to read, if there are tracers.
    - rLayer (int): The layer number to read, otherwise all layers are read.
    - finfo (bool): Whether to print information about the file.

    Returns:
    - F (np.ndarray): The field.
    - IDM (int): The number of grid points in the x-direction.
    - JDM (int): The number of grid points in the y-direction.
    - ll (int): The number of layers.
    '''
    try: 
        fgb = open(finb,'r')
    except:
        print('Could not open '+ finb)

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
        sys.exit('No idm found: Reading ', finb)

    IDM = int(data[0])
    data = fgb.readline().split()
    JDM = int(data[0])
    IJDM = IDM*JDM

    npad =4096-IJDM%4096

    if finfo:
        print('Reading HYCOM :{0} '.format(finb))

    aa = fgb.readline().split()

    cntr= 0
    nf = len(fld)
    FLOC=[]
    while cntr<1e6:
        aa = fgb.readline().split()
        if len(aa) == 0:
            break
        cntr += 1
        aname = aa[0]
        ii = aname.find(fld)
        if ii >= 0 and len(aname)==nf:
            FLOC.append(cntr)

    fgb.close()

    nrec = len(FLOC)
    if nrec == 0:
        raise Exception('read_ab: Field {0} not found in {1}'.format(fld, finb))

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
        nTR = dindx[0]+1
        nVlev = ll/nTR

    if nTR != 1 and finfo:
        print('read_ab: Found {0} variables {1}  per layer'.format(nTR, fld))

    lr1=-1
    lr2=-1
    if Rtrc is not None:
        if nTR < Rtrc:
            raise Exception('Number of saved tracers {0} < requested {1}'.\
                        format(nTR, Rtrc))
        dmm = np.copy(FLOC)
        FLOC = dmm[Rtrc-1::nTR]

        if lr1 < 0 or lr2 < 0 :
            lr2 = FLOC.shape[0]
            ll = lr2

    if rLayer is not None:
        lr1 = rLayer
        lr2 = lr1

    if lr1 < 0 or lr2 < 0:
        lr1 = 1
        lr2 = ll

    fga = open(fina,'rb')
    F = []
    ccL = -1
    huge = 0.0001*2.**99
    for ii in range(lr1,lr2+1):
        fga.seek(0)
        k0 = FLOC[ii-1]-1
        seek_position = np.int64(k0) * (npad + IJDM) * 4 
        fga.seek(seek_position, 0)
        dmm = np.fromfile(fga, dtype='>f',count=IJDM)
        amin = np.min(dmm[np.where(dmm<huge)])
        amax = np.max(dmm[np.where(dmm<huge)])
        if finfo:
            print('Reading {0} k={1} min={2:12.4f} max={3:12.4f}'.\
                format(fld,ii,amin,amax))
        dmm = dmm.reshape((JDM,IDM))
        ccL += 1

        if ccL == 0:
            F = np.copy(dmm)
        else:
            F = np.dstack((F,dmm))

    if ll == 0:
        print('!!! read_ab: {0} not found in {1} ERR'.format(fld, fina))
        print('!!! read hycom: check fields in ', finb)

    fga.close()

    return F, IDM, JDM, ll
######################################################
def find_indx_lonlat(x0, y0, LON, LAT):
    """
    Find closest grid point to lon/lat coordinate
    x0, y0 - geographic coordinates
    """
    if x0 > 180.:
        x0 = x0-360.

    dmm = np.sqrt((LON-x0)**2+(LAT-y0)**2)
    jj0, ii0 = np.where(dmm == np.min(dmm))

    return ii0[0], jj0[0]

def dx_dy(LON, LAT):
    """
    Find horizontal grid spacing from LON,LAT 2D arrays 
    hycom grid
    """
    IDM = LON.shape[1]
    JDM = LON.shape[0]
    print('Calculating DX, DY for HYCOM idm={0}, jdm={1}'.format(IDM,JDM))
    
    DX = np.zeros((JDM,IDM))
    DY = np.zeros((JDM,IDM))

    for ii in range(IDM-1):
        LT1 = LAT[:,ii]
        LT2 = LAT[:,ii+1]
        LN1 = LON[:,ii]
        LN2 = LON[:,ii+1]
        dx  = dist_sphcrd(LT1,LN1,LT2,LN2)
        DX[:,ii] = dx

    DX[:,IDM-1] = dx

    for jj in range(JDM-1):
        LT1 = LAT[jj,:]
        LT2 = LAT[jj+1,:]
        LN1 = LON[jj,:]
        LN2 = LON[jj+1,:]
        dy  = dist_sphcrd(LT1,LN1,LT2,LN2)
        DY[jj,:] = dy
        
    DY[JDM-1,:] = dy

    print('Min/max DX, m = {0:12.4f} / {1:12.4f}'.format(np.min(DX),np.max(DX)))
    print('Min/max DY, m = {0:12.4f} / {1:12.4f}'.format(np.min(DY),np.max(DY)))

    return DX, DY

def get_zz_zm(fina, finb, f_btm=True):
    """
    Derive layer depths: ZZ - interface depths,
    ZM - mid cell depths
    f_btm = true - ZZ=nan below bottom
    """
    rg = 9806.
    print(' Deriving ZZ, ZM from '+ fina)

    fld = 'thknss'
    F, nn, mm, ll = read_ab(fina,finb,fld,rLayer=1)
    F = F/rg
    F[np.where(F>1.e20)] = np.nan
    F[np.where(F<0.001)] = 0.
    ZZ = np.zeros((ll+1,mm,nn))
    ZZ[1,:,:] = -F
    ZM = np.zeros((ll,mm,nn))

    for kk in range(1,ll):
        F, nn, mm, ll = read_ab(fina, finb, fld, rLayer=kk+1)
        F = F/rg
        F[np.where(F>1.e20)] = np.nan
        F[np.where(F<0.001)] = 0.
        ZZ[kk+1,:,:] = ZZ[kk,:,:]-F

        if f_btm:
            jb, ib = np.where(F<0.001)
            ZZ[kk+1,jb,ib] = np.nan
    
    for kk in range(ll):
        ZM[kk,:,:] = 0.5*(ZZ[kk+1,:,:]+ZZ[kk,:,:])

    return ZZ, ZM

#####################################################################

import numpy as np
import matplotlib.pyplot as plt

# 1. Get u-vel data from the .a and .b files
U, IDM, JDM, _ = read_ab(a, b, 'u-vel.')

# 2. then find all of the gridpoints within a bounding box defined by the lon/lat coordinates (-90, 17) and (-82, 24)
lon_min, lon_max = -90, -82
lat_min, lat_max = 17, 24



# 3. then derive the layer depths for those gridpoints


# 4. then calculate the horizontal grid spacing for those gridpoints


# 5. then plot the u velocities for the first layer of the gridpoints


