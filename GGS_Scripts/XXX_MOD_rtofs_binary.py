# =========================
# X - Imports
# =========================

import numpy as np
import sys

# =========================
# RTOFS/HYCOM BINARY (DEPRECATED)
# =========================

### CLASS:
class RTOFS_binary:

    ### FUNCTION:
    def __init__(self, fina, finb, fld, Rtrc=None, rLayer=None, finfo=True):
        self.fina = fina
        self.finb = finb
        self.fld = fld
        self.Rtrc = Rtrc
        self.rLayer = rLayer
        self.finfo = finfo
        self.IDM = None
        self.JDM = None
        self.ll = None

    ### FUNCTION:
    def read_ab(self):
        try: 
            fgb = open(self.finb,'r')
        except:
            print('Could not open '+ self.finb)

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
            sys.exit('No idm found: Reading ', self.finb)

        IDM = int(data[0])
        data = fgb.readline().split()
        JDM = int(data[0])
        IJDM = IDM*JDM

        npad =4096-IJDM%4096

        if self.finfo:
            print('Reading HYCOM :{0} '.format(self.finb))

        aa = fgb.readline().split()

        cntr= 0
        nf = len(self.fld)
        FLOC=[]
        while cntr<1e6:
            aa = fgb.readline().split()
            if len(aa) == 0:
                break
            cntr += 1
            aname = aa[0]
            ii = aname.find(self.fld)
            if ii >= 0 and len(aname)==nf:
                FLOC.append(cntr)

        fgb.close()

        nrec = len(FLOC)
        if nrec == 0:
            raise Exception('read_ab: Field {0} not found in {1}'.format(self.fld, self.finb))

        ll = len(FLOC)
        if self.finfo:
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

        if nTR != 1 and self.finfo:
            print('read_ab: Found {0} variables {1}  per layer'.format(nTR, self.fld))

        lr1=-1
        lr2=-1
        if self.Rtrc is not None:
            if nTR < self.Rtrc:
                raise Exception('Number of saved tracers {0} < requested {1}'.\
                            format(nTR, self.Rtrc))
            dmm = np.copy(FLOC)
            FLOC = dmm[self.Rtrc-1::nTR]

            if lr1 < 0 or lr2 < 0 :
                lr2 = FLOC.shape[0]
                ll = lr2

        if self.rLayer is not None:
            lr1 = self.rLayer
            lr2 = lr1

        if lr1 < 0 or lr2 < 0:
            lr1 = 1
            lr2 = ll

        fga = open(self.fina,'rb')
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
            if self.finfo:
                print('Reading {0} k={1} min={2:12.4f} max={3:12.4f}'.\
                    format(self.fld,ii,amin,amax))
            dmm = dmm.reshape((JDM,IDM))
            ccL += 1

            if ccL == 0:
                F = np.copy(dmm)
            else:
                F = np.dstack((F,dmm))

        if ll == 0:
            print('!!! read_ab: {0} not found in {1} ERR'.format(self.fld, self.fina))
            print('!!! read hycom: check fields in ', self.finb)

        fga.close()

        return F, IDM, JDM, ll

    ### FUNCTION:
    def find_indx_lonlat(self, x0, y0, LON, LAT):
        if x0 > 180.:
            x0 = x0-360.

        dmm = np.sqrt((LON-x0)**2+(LAT-y0)**2)
        jj0, ii0 = np.where(dmm == np.min(dmm))

        return ii0[0], jj0[0]
    
    ### FUNCTION:
    def dx_dy(self, LON, LAT):
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
    
    ### FUNCTION:
    def get_zz_zm(self, fina, finb, f_btm=True):
        rg = 9806.
        print(' Deriving ZZ, ZM from '+fina)

        fld = 'thknss'
        F, nn, mm, ll = self.read_ab(fina,finb,fld,rLayer=1)
        F = F/rg
        F[np.where(F>1.e20)] = np.nan
        F[np.where(F<0.001)] = 0.
        ZZ = np.zeros((ll+1,mm,nn))
        ZZ[1,:,:] = -F
        ZM = np.zeros((ll,mm,nn))

        for kk in range(1,ll):
            F, nn, mm, ll = self.read_ab(fina,finb,fld,rLayer=kk+1)
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

### FUNCTION:
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

### FUNCTION:
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

### FUNCTION:
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

### FUNCTION:
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

# =========================
# RTOFS/HYCOM - Calculation Functions
# =========================

### FUNCTION:
def dist_sphcrd(xla1, xlo1, xla2, xlo2, Req=6371.0e3, Rpl=6357.e3):

    '''
    Calculates the great-circle distance between two geographical locations on an ellipse using Lambert formula

    Args:
    - xla1 (float): first point coordinates (latitude, longitude)
    - xlo1 (float): first point coordinates (latitude, longitude)
    - xla2 (float): second point coordinates (latitude, longitude)
    - xlo2 (float): second point coordinates (latitude, longitude)
    - Req (float): radius of the earth
    - Rpl (float): radius of the earth

    Returns:
    - dist_lmbrt (float): distance (in m)
    '''

    xla1 = np.float64(xla1)
    xlo1 = np.float64(xlo1)
    xla2 = np.float64(xla2)
    xlo2 = np.float64(xlo2)

    if np.absolute(xla1).max() > 90.0:
        print("ERR: dist_sphcrd Lat1 > 90")
        dist = float("nan")
        return dist
    if np.absolute(xla2).max() > 90.0:
        print("ERR: dist_sphcrd Lat2 > 90")
        dist = float("nan")
        return dist

    cf = np.pi/180.
    phi1 = xla1*cf
    phi2 = xla2*cf
    lmb1 = xlo1*cf
    lmb2 = xlo2*cf
    dphi = abs(phi2-phi1)
    dlmb = abs(lmb2-lmb1)

    I0s = []

    if isinstance(dphi,float) or isinstance(dphi,int):
        if dphi == 0.0 and dlmb == 0.0:
            dist_lmbrt = 0.0
            return dist_lmbrt
    else:
        if np.min(dphi) == 0.0 and np.min(dlmb[dphi==0])==0:
            I0s = np.argwhere((dphi==0.0) & (dlmb==0.0))
            if len(I0s)>0:
                dphi[I0s]=1.e-8
                dlmb[I0s]=1.e-8

    Eflat = (Req-Rpl)/Req
    aa1 = (np.sin(dphi/2.))**2
    aa2 = np.cos(phi1)*np.cos(phi2)*(np.sin(dlmb/2.))**2
    dsgm_hv = 2.*np.arcsin(np.sqrt(aa1+aa2))

    beta1 = np.arctan((1.-Eflat)*np.tan(phi1))
    beta2 = np.arctan((1.-Eflat)*np.tan(phi2))
    PP = 0.5*(beta1+beta2)
    QQ = 0.5*(beta2-beta1)
    X = (dsgm_hv-np.sin(dsgm_hv))*( (np.sin(PP))**2*(np.cos(QQ))**2 )/( (np.cos(dsgm_hv/2.))**2 )
    Y = (dsgm_hv+np.sin(dsgm_hv))*( (np.cos(PP))**2*(np.sin(QQ))**2 )/( (np.sin(dsgm_hv/2.))**2 )

    dist_lmbrt = Req*(dsgm_hv-Eflat/2.*(X+Y))

    if np.min(dist_lmbrt)<0.0:
        print('WARNING: spheric distance <0: ', np.min(dist_lmbrt))

    return dist_lmbrt