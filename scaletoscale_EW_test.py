#!/usr/bin/env python
#################################################
# Load Modules
#################################################
from netCDF4 import Dataset
import numpy as np
from scipy.ndimage import filters as filt

from functools import partial

import sys, os
import niskindata
import R_tools_new_goth as tN
# import R_tools_fort as tF

def ncopen(dr, filestr, depth, opt = 'r'):
    return Dataset(dr + filestr+'.{0:04}'.format(depth)+'.nc', opt)

######################################

def makenan(x, val = np.nan):
    x[x==0] = val 
    x[x<-100000] = val 
    x[x>1000000] = val

######################################
# modified by zsm
# def loadvd(nch, ncz, itime):
#     vd = {}
#     for var in ['u', 'v', 'w', 'uz', 'vz']:
#         if var == 'u':
#             func = tN.u2rho
#         elif var == 'v':
#             func = tN.v2rho
#         else:
#             func = None
#         try:
#                 vd[var] = np.squeeze(tN.ncload(nch, var, itime, func = func))
#         except: 
#                 vd[var] = np.squeeze(tN.ncload(ncz, var, itime, func = func))
#         makenan(vd[var]) 
#         #print(vd[var].shape)
#     return vd

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
                vd[var] = np.squeeze(tN.ncload(nch, var, itime, func = func))
        except: 
                vd[var] = np.squeeze(tN.ncload(ncz, var, itime, func = func))
        makenan(vd[var]) 
        #print(vd[var].shape)
    return vd
######################################
# modified by zsm
# def loadvdEW(nch, itime):
#     vde = {}
#     vdw = {}
#     for var in ['u', 'v', 'w', 'uz', 'vz']:
#         vde[var] = np.squeeze(tN.ncload(nch, var  + 'Eddy', itime))
#         vdw[var] = np.squeeze(tN.ncload(nch, var  + 'Wave', itime))
#         #makenan(vd[var]) 
#         print(vd[var].shape)
#     return vde, vdw

def loadvdEW(nch, itime):
    vde = {}
    vdw = {}
    for var in ['u', 'v']:
        vde[var] = np.squeeze(tN.ncload(nch, var  + 'Eddy', itime))
        vdw[var] = np.squeeze(tN.ncload(nch, var  + 'Wave', itime))
        #makenan(vd[var]) 
        print(vd[var].shape)
    return vde, vdw
######################################
def filtr(var, n):
    return filt.uniform_filter(var,n,None,'mirror',0)
        
######################################
def prodf(u,v,n):
    return filtr(u*v,n)-filtr(u,n)*filtr(v,n)

######################################
def tau(dims,vd,n,vd2 = None):
    #print conv[dims[0]], conv[dims[1]]
    if vd2 is None:
        return prodf(vd[dims[0]], vd[dims[1]], n)
    else:
        return [prodf(vd[dims[0]], vd2[dims[1]], n), prodf(vd2[dims[0]], vd[dims[1]], n)]
        
######################################
def tau2(dims,vd,n,vd2 = None):
    #print conv[dims[0]], conv[dims[1]]
    if vd2 is None:
        return prodf(vd[dims[0]], vd[dims[1]], n)
    else:
        return prodf(vd[dims[0]], vd2[dims[1]], n)
######################################
        
def mapfilt(varlist, n):
    return [filtr(var,n) for var in varlist]

# Roy modified
#Narray = list(range(1,12))
#Narray.extend(list(range(12,int(110/scalefac),int(6/scalefac))))
#Narray = np.array(Narray)[0:-1:2]

#def nanmeanzones(varval, Nbands):
    
#w->div, e->rot
def tauijVij(vd, vdw, vde, gd, varlist):
    r = tN.Zeros(N0)
    #Tm = tN.Zeros(N0)
    od = {}
    for var in varlist:
        od[var] = tN.Zeros(N0)
    #print('Starting flux loop with %d steps' %(N0))
    #print('Filter order:', end=' ')
    print('Starting flux loop with %d steps' %(N0))
    print('Filter order:', end=' ')
    for I, n in enumerate(Narray):
        print('%d ' %n, end=' ')
        r[I] = n*gd['dx']
        
    #w->div, e->rot
        ux, uy, vx, vy = mapfilt(list(tN.shear2d(vd['u'],vd['v'],gd['pm'],gd['pn'])), n)
        uxe, uye, vxe, vye = mapfilt(list(tN.shear2d(vde['u'],vde['v'],gd['pm'],gd['pn'])), n)  
        uxw, uyw, vxw, vyw = mapfilt(list(tN.shear2d(vdw['u'],vdw['v'],gd['pm'],gd['pn'])), n)  
        
        
        #print 'Done smoothing'
        tauxx, tauyy, tauxy, tauyx = [tau('uu', vd, n), tau('vv', vd, n), tau('uv', vd, n), tau('vu', vd, n)]   
        tauxxww, tauyyww, tauxyww, tauyxww = [tau('uu', vdw, n), tau('vv', vdw, n), tau('uv', vdw, n), tau('vu', vdw, n)]   
        tauxxee, tauyyee, tauxyee, tauyxee = [tau('uu', vde, n), tau('vv', vde, n), tau('uv', vde, n), tau('vu', vde, n)]   
        tauxxew, tauyyew, tauxyew, tauyxew = [tau2('uu', vde, n, vdw), tau2('vv', vde, n,vdw), tau2('uv', vde, n, vdw), tau2('vu', vde, n, vdw)]    
        tauxxwe, tauyywe, tauxywe, tauyxwe = [tau2('uu', vdw, n, vde), tau2('vv', vdw, n,vde), tau2('uv', vdw, n, vde), tau2('vu', vdw, n, vde)]    

        SPH    = -(tauxx*ux    + tauxy*uy    + tauyx*vx    + tauyy*vy) 
        SPHeeE = -(tauxxee*uxe + tauxyee*uye + tauyxee*vxe + tauyyee*vye) 
        SPHwwW = -(tauxxww*uxw + tauxyww*uyw + tauyxww*vxw + tauyyww*vyw) 
        SPHeeW = -(tauxxee*uxw + tauxyee*uyw + tauyxee*vxw + tauyyee*vyw) 
        SPHwwE = -(tauxxww*uxe + tauxyww*uye + tauyxww*vxe + tauyyww*vye) 
        SPHweE = -(tauxxwe*uxe + tauxywe*uye + tauyxwe*vxe + tauyywe*vye) 
        SPHewE = -(tauxxew*uxe + tauxyew*uye + tauyxew*vxe + tauyyew*vye) 
        SPHweW = -(tauxxwe*uxw + tauxywe*uyw + tauyxwe*vxw + tauyywe*vyw) 
        SPHewW = -(tauxxew*uxw + tauxyew*uyw + tauyxew*vxw + tauyyew*vyw) 

        # modifiedy by zsm, calc vertical term
        # uze, vze, uzw, vzw, uz, vz = mapfilt([vde['uz'], vde['vz'], vdw['uz'], vdw['vz'], vd['uz'], vd['vz']], n)    
        # taux, tauy = [tau('uw', vd, n), tau('vw', vd, n)] 
        # tauxww, tauyww = [tau('uw', vdw, n), tau('vw', vdw, n)] 
        # tauxee, tauyee = [tau('uw', vde, n), tau('vw', vde, n)] 
        # tauxew, tauxwe = tau('uw', vde, n, vdw)
        # tauyew, tauywe = tau('vw', vde, n, vdw)
        ## tauxewe = tauxew +  tauxwe
        ## tauyewe = tauyew +  tauywe

        # SPVeeE = -(tauxee*uze + tauyee*vze)
        # SPVwwW = -(tauxww*uzw + tauyww*vzw)
        # SPVeeW = -(tauxee*uzw + tauyee*vzw)
        # SPVwwE = -(tauxww*uze + tauyww*vze)
        ##SPVwewW = -(tauxewe*uzw + tauyewe*vzw)
        ##SPVwewE = -(tauxewe*uze + tauyewe*vze)   
        # SPVewE = -(tauxew*uze + tauyew*vze)   
        # SPVweE = -(tauxwe*uze + tauywe*vze)   
        
        # SPVewW  = -(tauxew*uzw + tauyew*vzw)
        # SPVweW  = -(tauxwe*uzw + tauywe*vzw)
        # SPV  = -(taux  *uz  + tauy  *vz)    
        
        # SPVtot = SPVeeE +SPVwwW +SPVeeW +SPVwwE +SPVewW + SPVewE + SPVweW + SPVweE
        SPHtot = SPHeeE + SPHwwW + SPHeeW + SPHwwE + SPHweE + SPHewE + SPHweW + SPHewW
        for var in varlist:
            od[var][I] = np.squeeze(np.nanmean(locals()[var]))
    return r, od 
######################################
#############MAIN#####################
######################################
#Output
# varlist = [ 'SPHeeE', 'SPHeeW', 'SPHwwW', 'SPHwwE', 'SPHweW', 'SPHewW', 'SPHewE', 'SPHweE',
#             'SPVeeE', 'SPVeeW', 'SPVwwW', 'SPVwwE', 'SPVweW', 'SPVewW', 'SPVewE', 'SPVweE',
#            'SPH', 'SPHtot', 'SPV', 'SPVtot']
varlist = [ 'SPHeeE', 'SPHeeW', 'SPHwwW', 'SPHwwE', 'SPHweW', 'SPHewW', 'SPHewE', 'SPHweE',
           'SPH', 'SPHtot']
####################################################################################
try:
        grd = sys.argv[1] 
        season = sys.argv[2]
        depth = int(sys.argv[3])
        simul = sys.argv[4]
        opt = sys.argv[5]
except:
        print('Usage python scaletoscale_EW.py grd[2km/500] season depth simul[hf/smooth] opt[helm/nohelm]')
        sys.exit(0) 
dr, dro, gridfile, filehis, filez = niskindata.path(grd, season, simul)
gd = tN.gridDict(dr,gridfile)
nch = ncopen(dr, filehis, depth)
# modified by zsm
# ncz = ncopen(dr, filez, depth)
ncew = ncopen(dro, 'eddyWave', depth)

filtE = ncew.filterEddyTime
filtW = ncew.filterWaveTime
outfile = 's2sEddyWave' if opt =='nohelm' else 's2sEddyWaveRotDiv'
outfile += str(filtE) + '_' + str(filtW)

nco = ncopen(dro, 's2sEddyWave', depth, 'w')
Ny = gd['Ny']
print(Ny)

Nt = nch.dimensions['time'].size
Nt0 = 0
######################################
if grd == '2km': 
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
print('Narray:', np.array(Narray) * gd['dx']/1000, len(Narray))
######################################
nco.createDimension('filtscale', N0)
nco.createDimension('time', None)
for var in varlist: 
    nco.createVariable(var,np.dtype('float32').char,('time', 'filtscale'))
nco.createVariable('filtscale',np.dtype('float32').char,('filtscale'))

atime = 0
for itime in range(Nt0,Nt,1):
    print(f'itime: {itime}')
    # modified by zsm
    # vd  = loadvd(nch, ncz, itime)
    vd  = loadvd(nch, nch, itime)
    #w->div, e->rot
    vde, vdw = loadvdEW(ncew, itime)
    filtscale, od = tauijVij(vd, vdw, vde, gd, varlist)
    for var in varlist: 
        nco.variables[var][atime, :] = od[var][:]
    if itime==0:
        nco.variables['filtscale'][:] = filtscale
    atime+=1    
nco.close()
nch.close()
# nchelm.close()
# ncz.close()
