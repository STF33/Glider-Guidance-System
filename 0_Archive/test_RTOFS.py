# Imports
import numpy as np
import sys
import numpy as np
import os
from netCDF4 import Dataset
from scipy.interpolate import griddata

# Paths to the .a and .b files
a = 'C:/Users/salfr/Downloads/GGS_TEST/rtofs_glo.t00z.f12.archv.a'
b = 'C:/Users/salfr/Downloads/GGS_TEST/rtofs_glo.t00z.f12.archv.b'

# Function to read the .a and .b files
def read_rtofs(fina, finb, fld, Rtrc=None, rLayer=None, finfo=True):

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

    npad =4096-IJDM%4096

    if finfo:
        print('Reading HYCOM :{0} '.format(finb))

    aa = fgb.readline().split()

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
        raise Exception('read_rtofs: Field {0} not found in {1}'.format(fld,finb))

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
        nTR = dindx[0]+1
        nVlev = ll/nTR

    if nTR != 1 and finfo:
        print('read_rtofs: Found {0} variables {1}  per layer'.format(nTR,fld))

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
        dmm = np.fromfile(fga, dtype='>f',count=IJDM) # read 1 layer
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
        print('!!! read_rtofs: {0} not found in {1} ERR'.format(fld,fina))
        print('!!! read hycom: check fields in ',finb)

    fga.close()

    return F, IDM, JDM, ll

# Run the function: 'read_rtofs'
F, IDM, JDM, ll = read_rtofs(a, b, 'thknss')
F, IDM, JDM, ll = read_rtofs(a, b, 'u-vel.')
F, IDM, JDM, ll = read_rtofs(a, b, 'v-vel.')


