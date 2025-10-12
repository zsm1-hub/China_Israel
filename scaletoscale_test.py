#!/usr/bin/env python
#################################################
# Load Modules
#################################################

from netCDF4 import Dataset
import numpy as np
from scipy.ndimage import filters as filt

from functools import partial

import sys, os

import R_tools_new_goth as tN
import niskindata

print('Loading netcdf files\n')
def ncopen(dr, filestr, depth, opt = 'r'):
    filen = dr + filestr+'.{0:04}'.format(depth)+'.nc'
    print('*******************************************')
    print(filen)
    return Dataset(filen, opt)

grd = sys.argv[1]
season = sys.argv[2]
depth = int(sys.argv[3])
Case = (sys.argv[4])
dr, _, gridfile, filehis, filez = niskindata.path(grd, season, Case)
                
######################################
# MAIN
#####################################
opt = 'mirror'
######################################


gd = tN.gridDict(dr,gridfile)
# ncz = ncopen(dr, filez, depth)

nch = ncopen(dr, filehis, depth)
# nchelm = ncopen(dr, 'helmholtz', depth)
nco = ncopen(dr, 's2sflux_spec', depth, 'w')
ncoz = ncopen(dr, 's2sflux_spec_all_'+Case, depth, 'w')

Nx = nch.dimensions['xi_rho' ].size
Ny = nch.dimensions['eta_rho' ].size
assert(Nx==Ny)

####################################

def makenan(x, val = np.nan):
    x[x==0] = val 
    x[x<-100000] = val 
    x[x>1000000] = val

######################################
def loadvd_h(nch, itime):
    vd = {}
    vd['u'] = np.squeeze(tN.ncload(nch, 'udiv', itime, func = tN.u2rho)) 
    vd['v'] = np.squeeze(tN.ncload(nch, 'vdiv', itime, func = tN.v2rho)) 
    return vd

def loadvd(nch, ncz, itime):
    vd = {}
    for var in ['u', 'v']:
        if var == 'u':
            func = tN.u2rho
        elif var == 'v':
            func = tN.v2rho
        else:
            func = None
        try:
                varval = np.squeeze(tN.ncload(nch, var, itime, func = func))
                vd[var] = varval[:,:]
        except: 
                varval = np.squeeze(tN.ncload(ncz, var, itime, func = func))
                vd[var] = varval[:,:]
        makenan(vd[var]) 
        print(vd[var].shape)
    return vd
######################################
from scipy.fft import rfft2, irfft2
from skimage.filters import window
from scipy.signal import get_window


kx = np.fft.fftfreq(Nx, d=1./Nx)
ky = np.fft.fftfreq(Ny, d=1./Ny)
nk = Nx//2 + 1
KY, KX = np.meshgrid(ky[:nk], kx)

def sfilt(ur, r):
    # w = window(('tukey', 0.1), ur.shape)
    # u = rfft2(ur*w)
    u = rfft2(ur)
    nx, nk = u.shape
    # modified by zsm, effective scale
    # nc = nx/2.0/r
    nc = nx/r
    K2 = KX**2 + KY**2
    mask = K2>nc**2
    u[mask] = 0
    return irfft2(u)

def sfilter_corr(ur, r, pad_ratio=0.1, bry='periodic'):
    """
    mirrors
    periodic
    """
    if bry=='mirrors':
        Mode='reflect'
    if bry=='periodic':
        Mode='wrap'
    Ny, Nx = ur.shape
    
    # 1. calc size of mirrors
    # pad_x = int(Nx * pad_ratio)
    # pad_y = int(Ny * pad_ratio)
    pad_x = int(r)
    pad_y = int(r)
    
    # 2. mirrors
    ur_padded = np.pad(ur, ((pad_y, pad_y), (pad_x, pad_x)), mode=Mode)
    Ny_pad, Nx_pad = ur_padded.shape
    
    # 3. pad
    if bry=='mirrors':
        window_x = get_window(('tukey', 0.1), Nx_pad)
        window_y = get_window(('tukey', 0.1), Ny_pad)
    # if bry=='periodic':
        # window_x = get_window(('hann'), Nx_pad)
        # window_y = get_window(('hann'), Ny_pad)
        # window_x = get_window(('tukey', 0.1), Nx_pad)
        # window_y = get_window(('tukey', 0.1), Ny_pad)
        # window_x = get_window(('kaiser', 14), Nx_pad)
        # window_y = get_window(('kaiser', 14), Ny_pad)
    #     window_x = np.ones((Nx_pad))
    #     window_y = np.ones((Ny_pad))
    
    # w = np.outer(window_y, window_x)
    
    # 4. add window fft
    # u = rfft2(ur_padded*w)
    u = rfft2(ur_padded)
    
    
    # 5. 
    kx_pad = np.fft.fftfreq(Nx_pad, d=1./Nx_pad)  # 
    ky_pad = np.fft.fftfreq(Ny_pad, d=1./Ny_pad)
    nk_pad = Nx_pad // 2 + 1
    
    # build grid
    KX_pad, KY_pad = np.meshgrid(kx_pad[:nk_pad], ky_pad)
    
    # 6. calc cutoff wavenumber
    nc = Nx / (2.0 * r)  # 
    
    # 7. mask
    K2 = KX_pad**2 + KY_pad**2
    mask = K2 > nc**2
    
    # 8. low pass
    u[mask] = 0
    
    # 9. ifft
    ur_filtered_padded = irfft2(u, s=(Ny_pad, Nx_pad))
    ur_filtered = ur_filtered_padded[pad_y:pad_y+Ny, pad_x:pad_x+Nx]
    
    return ur_filtered
 
def filtr(var, n):
    return sfilt(var, n)
# def filtr(var, n):
#    return filt.uniform_filter(var,n,None,opt,0)
# def filtr(var, n):
#     return sfilter_corr(var, n)        
######################################
def prodf(u,v,n):
    return filtr(u*v,n)-filtr(u,n)*filtr(v,n)

######################################
def tau(dims,vd,n):
    conv = {'x':'u', 'y':'v', 'z':'w'}
    return prodf(vd[conv[dims[0]]], vd[conv[dims[1]]], n)

######################################
def mapfilt(varlist, n):
    return [filtr(var,n) for var in varlist]

def jac(ux, uy, vx, vy):
    return ux * vy - vx * uy

if grd == '2km'or grd =='2km_smooth': 
    Narray = list(range(1,18))
    Narray.extend(list(range(18,int(48),int(3))))
    Narray.extend(list(range(48,int(120),int(6))))
    N0 = len(Narray)
else:
    Narray = list(range(1,18))
    Narray.extend(list(range(18,int(48),int(3))))
    Narray.extend(list(range(48,int(120),int(6))))
    Narray.extend(list(range(120,int(240),int(12))))
    Narray.extend(list(range(240,int(360),int(24))))
    N0 = len(Narray)

if Case=='hit':
    Narray = list(range(1,18))
    Narray.extend(list(range(18,int(48),int(3))))
    Narray.extend(list(range(48,int(120),int(6))))
    Narray.extend(list(range(120,int(240),int(12))))
    Narray.extend(list(range(240,int(520),int(24))))
    N0 = len(Narray)

if Case=='aviso':
    Narray = list(range(1,18))
    Narray.extend(list(range(18,int(48),int(3))))
    Narray.extend(list(range(48,int(120),int(6))))
    Narray.extend(list(range(120,int(240),int(12))))
    N0 = len(Narray)


print('Narray:', np.array(Narray) * gd['dx']/1000, len(Narray))


######################################
def filtflux(vd, gd):
    r = tN.Zeros(N0)
    Tm = tN.Zeros(N0)
    Thm = tN.Zeros(N0)
    Tzm = tN.Zeros(N0)
    Tdeltam = tN.Zeros(N0)
    Talpham = tN.Zeros(N0)
    Pm = tN.Zeros(N0)
    Ekm = tN.Zeros(N0)
    deltam = tN.Zeros(N0)    
    alpham = tN.Zeros(N0)
    alpharm = tN.Zeros(N0)
    alphadm = tN.Zeros(N0)
    alphardm = tN.Zeros(N0)
    Th_zsm=tN.Zeros([Nx,Ny,N0])
    print('Starting flux loop with %d steps' %(N0))
    print('Filter order:', end=' ')
    for I, n in enumerate(Narray):
        #Filter scale, number of points in filter is n*Nf
        print('%d ' %n, end=' ')
        r[I] = n*gd['dx']
        ux, uy, vx, vy = tN.shear2d(vd['u'],vd['v'],gd['pm'],gd['pn'])  
        # uxd, uyd, vxd, vyd = tN.shear2d(vdd['u'],vdd['v'],gd['pm'],gd['pn'])  
        
        #print 'Done shear'
        uxr, uyr, vxr, vyr = mapfilt([ux, uy, vx, vy], n)    
        # uxd, uyd, vxd, vyd = mapfilt([uxd, uyd, vxd, vyd], n)

        #print 'Done smoothing'
        # delta = vyr + uxr
        # zeta = vxr - uyr
        
        # alpha = (uxr - vyr)**2 + (vxr + uyr)**2
        # alphar = zeta**2 - 4*jac(uxr - uxd, uyr - uyd, vxr - vxd, vyr - vyd)
        # alphad = delta**2 - 4*jac(uxd, uyd, vxd, vyd)
        # alphard = -(4*jac(uxd, uyr - uyd, vxd, vyr - vyd)) 
        # alphard -= (4*jac(uxr - uxd, uyd, vxr - vxd, vyd))

        Th = -(tau('xy',vd, n) * (uyr + vxr)  +  tau('xx', vd, n) * uxr + tau('yy', vd, n) * vyr )
        # Tz = -(tau('xz', vd, n) * uzr + tau('yz', vd, n) * vzr)
        # P = prodf(vd['w'],vd['b'],n)
        Ek = 0.5*(tau('xx', vd, n) + tau('yy', vd, n))
        # Tdelta = -Ek*delta
        #print 'Done with T and P '
        Thm[I], Ekm[I]  = list(map(np.nanmean, [Th[:], Ek[:]]))
        # Th_zsm[:,:,I]=tN.v2rho(Th)
        Th_zsm[:,:,I]=Th
        print('*******************zsm******************')
        print(Th_zsm.shape)
        print(uy.shape)
        # print(Narray)
        # alpham[I], alpharm[I], alphadm[I], alphardm[I], deltam[I] = list(map(np.nanmean,[alpha, alphar, alphad, alphard, delta]))
        
        # Tm[I] = Thm[I]+Tzm[I]
        # Talpham[I] = Thm[I] - Tdeltam[I]
    #print r, Tsm, P
    return r, Thm, Ekm, Th_zsm

#main
try:
        Nt = int(sys.argv[5])
except:
        Nt = nch.dimensions['time'].size
Nt0 = 0
print('eventually **********************************')
print(Nt)
#Nt = 1000 
varlist = ['Th', 'Ek']

varlist2 = []  # 创建空列表
for xsc in range(N0):
    # 使用 append() 添加元素
    varlist2.append('Th'+str(Narray[xsc]))
    
print(varlist2)

#Size of the maximum filter scale is Nf-th of the domain
atime = 0
#Output
for var in  varlist:
    globals()[var+'m'] = tN.Zeros(N0)

# file1
filtscale = tN.Zeros(N0)
nco.createDimension('filtscale', len(filtscale))
for var in varlist: 
    nco.createVariable(var+'m',np.dtype('float32').char,('filtscale'))
nco.createVariable('filtscale',np.dtype('float32').char,('filtscale'))
# file2
ncoz.createDimension('time', Nt)        # 无限制大小，时间维度
ncoz.createDimension('depth', 1)   # 深度维度
ncoz.createDimension('eta_rho', Ny)    # eta_rho 维度
ncoz.createDimension('xi_rho', Nx) 
for var in varlist2:
    ncoz.createVariable(var, 'f4', ('time', 'depth', 'eta_rho', 'xi_rho'))


for itime in range(Nt0,Nt,1):
    print('time:%d of %d' %(itime, Nt))
    
    #vd has the flow variables, gd the grid dimensions
    vd = loadvd(nch, nch, itime)
    # vdd = loadvd_h(nchelm, itime) 
    #r is the filter scale, the rest are fluxes averaged over (x,y) space and a function of r
    filtscale, Th, Ek,Th_zsm = filtflux(vd, gd)
    print(Th.shape, Ek.shape)
    for var in varlist:
        globals()[var+'m'] += globals()[var]/float(Nt-Nt0)               
    #create output file structures in first iteration
    nn=0
    for var in varlist2:
        ncoz.variables[var][itime,:,:,:]=Th_zsm[:,:,nn].T
        nn=nn+1
for var in varlist:
        nco.variables[var+'m'][:] = globals()[var+'m']
nco.variables['filtscale'][:] = filtscale
#print(var, globals()[var].shape, nco.variables[var])
    
nco.close()
ncoz.close()

#print 'Flux for r=%f km, T(r)=%e cm2s-3' %(r[I]/1000,10000*Tm[I])
