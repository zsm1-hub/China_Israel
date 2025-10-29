#!/usr/bin/env python
#################################################
# Load Modules
#################################################

from netCDF4 import Dataset
import numpy as np
from scipy.ndimage import filters as filt
import ray  # 导入 Ray

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

# 初始化 Ray
ray.init(num_cpus=39)

gd = tN.gridDict(dr,gridfile)
# ncz = ncopen(dr, filez, depth)

nch = ncopen(dr, filehis, depth)
# nchelm = ncopen(dr, 'helmholtz', depth)
nco = ncopen(dr, 's2sflux_spec_'+Case, depth, 'w')
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
    # w = window('hann', ur.shape)
    # Nx_pad, Ny_pad = ur.shape
    # window_x = get_window(('kaiser', 1), Nx_pad)
    # window_y = get_window(('kaiser', 1), Ny_pad)
    # w = np.outer(window_y, window_x)
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

def sfilt_corr(ur, r):
    # w = window(('tukey', 0.1), ur.shape)
    # w = window('hann', ur.shape)
    Nx_pad, Ny_pad = ur.shape
    window_x = get_window(('kaiser', 1), Nx_pad)
    window_y = get_window(('kaiser', 1), Ny_pad)
    w = np.outer(window_y, window_x)
    u = rfft2(ur*w)
    # u = rfft2(ur)
    nx, nk = u.shape
    # modified by zsm, effective scale
    # nc = nx/2.0/r
    nc = nx/r
    K2 = KX**2 + KY**2
    mask = K2>nc**2
    u[mask] = 0
    
    window_energy = np.sum(w**2)
    correction_factor = np.sqrt(w.size / window_energy)
    ucorr=tN.v2rho(irfft2(u))*correction_factor
    # ucorr=(irfft2(u))*correction_factor
    # ucorr=(irfft2(u))
    return ucorr

# def filtr(var, n):
#     return sfilt(var, n)
def filtr(var, n):
    return sfilt_corr(var, n)
# def filtr(var, n):
#    return filt.uniform_filter(var,n,None,opt,0)
    
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
    Narray = list(range(2,18))
    Narray.extend(list(range(18,int(48),int(3))))
    Narray.extend(list(range(48,int(120),int(6))))
    Narray.extend(list(range(120,int(240),int(12))))
    Narray.extend(list(range(240,int(520),int(24))))
    Narray.append(512)
    N0 = len(Narray)

if Case=='aviso':
    Narray = list(range(1,18))
    Narray.extend(list(range(18,int(48),int(3))))
    Narray.extend(list(range(48,int(120),int(6))))
    Narray.extend(list(range(120,int(240),int(12))))
    N0 = len(Narray)


print('Narray:', np.array(Narray) * gd['dx']/1000, len(Narray))


######################################
# 将filtflux函数改为Ray远程函数
@ray.remote
def filtflux(vd, gd, Narray, N0, Nx, Ny):
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
        uxr, uyr, vxr, vyr = mapfilt([ux, uy, vx, vy], n)   
        print('check if filt dimension is right !')
        
        Th = -(tau('xy',vd, n) * (uyr + vxr)  +  tau('xx', vd, n) * uxr + tau('yy', vd, n) * vyr )
        # P = prodf(vd['w'],vd['b'],n)
        Ek = 0.5*(tau('xx', vd, n) + tau('yy', vd, n))
        # Tdelta = -Ek*delta
        #print 'Done with T and P '
        Thm[I], Ekm[I]  = list(map(np.nanmean, [Th[:], Ek[:]]))
        # Th_zsm[:,:,I]=tN.v2rho(Th)
        Th_zsm[:,:,I]=Th
        print('*******************zsm******************')
        print(Th_zsm.shape)
        print(Th.shape)
        print(ux.shape)
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

# 存储所有任务的结果
results = []

# 提交任务
for itime in range(Nt0, Nt, 1):
    print('Submitting task for time:%d of %d' %(itime, Nt))
    vd = loadvd(nch, nch, itime)
    # 提交远程任务
    results.append(filtflux.remote(vd, gd, Narray, N0, Nx, Ny))

# 收集结果
all_results = ray.get(results)

# 初始化用于平均的累加器
Thm_sum = np.zeros(N0)
Ekm_sum = np.zeros(N0)
Th_zsm_all = np.zeros((Nt, Nx, Ny, N0))  # 注意：这里假设Th_zsm是每个时间步的每个尺度的二维场

# 遍历每个时间步的结果
for i, (r, Th, Ek, Th_zsm) in enumerate(all_results):
    # 累加Thm和Ekm
    Thm_sum += Th
    Ekm_sum += Ek
    # 存储每个时间步的Th_zsm
    for j in range(N0):
        var_name = 'Th'+str(Narray[j])
        # 将Th_zsm[:,:,j]存入ncoz文件中对应的时间步和变量
        ncoz.variables[var_name][i, 0, :, :] = Th_zsm[:, :, j].T  # 注意转置，如果维度顺序需要调整

# 计算平均
Thm_avg = Thm_sum / (Nt - Nt0)
Ekm_avg = Ekm_sum / (Nt - Nt0)

# 将平均结果写入nco文件
nco.variables['Thm'][:] = Thm_avg
nco.variables['Ekm'][:] = Ekm_avg
nco.variables['filtscale'][:] = r  # 注意：这里r是最后一个时间步的r，但实际上每个时间步的r相同，所以取任意一个都可以

nco.close()
ncoz.close()

# 关闭Ray
ray.shutdown()