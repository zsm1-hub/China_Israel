#!/usr/bin/env python
# coding: utf-8

# # Parcels environment : base

# In[17]:


import parcels
import numpy as np
from datetime import timedelta
from glob import glob
import matplotlib.pyplot as plt
import xarray as xr
import os

def u2rho_2d (var_u):
    [Mp,L]=var_u.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((Mp,Lp))
    var_rho[:,1:L-1]=0.5*(var_u[:,0:Lm-1]+var_u[:,1:L-1])
    var_rho[:,0]=var_rho[:,1]
    var_rho[:,Lp-1]=var_rho[:,-2]
    return var_rho
    
def v2rho_2d (var_v):
    [M,Lp]=var_v.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((Mp,Lp))
    var_rho[1:M-1,:]=0.5*(var_v[0:Mm-1,:]+var_v[1:M-1,:])
    var_rho[0,:]=var_rho[1,:]
    var_rho[Mp-1,:]=var_rho[-2,:]
    return var_rho

def u2rho_3d (var_u):
    [N,Mp,L]=var_u.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((N,Mp,Lp))
    var_rho[:,:,1:-1]=0.5*(var_u[:,:,1:]+var_u[:,:,:-1])
    var_rho[:,:,0]=var_rho[:,:,1]
    var_rho[:,:,-1]=var_rho[:,:,-2]
    return var_rho
    
def v2rho_3d (var_v):
    [N,M,Lp]=var_v.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((N,Mp,Lp))
    var_rho[:,1:-1,:]=0.5*(var_v[:,1:,:]+var_v[:,:-1,:])
    var_rho[:,0,:]=var_rho[:,1,:]
    var_rho[:,-1,:]=var_rho[:,-2,:]
    return var_rho

def spheric_dist(lat1, lat2, lon1, lon2):
    """
    Compute the spherical distance between two points on Earth.
    
    Parameters:
    lat1, lat2 : array-like
        Latitude of the two points (in degrees).
    lon1, lon2 : array-like
        Longitude of the two points (in degrees).
    
    Returns:
    dist : array-like
        The spherical distance between the points (in meters).
    """
    
    # Earth radius in meters
    R = 6367442.76
    
    # Determine proper longitudinal shift
    l = np.abs(lon2 - lon1)
    l[l >= 180] = 360 - l[l >= 180]
    
    # Convert decimal degrees to radians
    deg2rad = np.pi / 180
    lat1 = lat1 * deg2rad
    lat2 = lat2 * deg2rad
    l = l * deg2rad
    
    # Compute the distances
    dist = R * np.arcsin(np.sqrt(((np.sin(l) * np.cos(lat2)) ** 2) + 
                                 ((np.sin(lat2) * np.cos(lat1)) - 
                                  (np.sin(lat1) * np.cos(lat2) * np.cos(l))) ** 2))
    
    return dist

def spheric_dist_one(lat1, lat2, lon1, lon2):
    """计算两个经纬度点之间的球面距离（标量版）"""
    # 处理经度差
    l = np.abs(lon2 - lon1)
    if l >= 180:
        l = 360 - l
    
    # 转换为弧度
    deg2rad = np.pi / 180
    lat1_rad = lat1 * deg2rad
    lat2_rad = lat2 * deg2rad
    l_rad = l * deg2rad
    
    # 球面距离公式
    distance = 6371 * np.arccos(
        np.sin(lat1_rad) * np.sin(lat2_rad) + 
        np.cos(lat1_rad) * np.cos(lat2_rad) * np.cos(l_rad)
    )
    return distance
def transunit_spher2flat(u,v,lat):
    v1=v*1852*60
    u1=u*1852*60*np.cos(lat*np.pi/180)
    return u1,v1

def trans_vel_roms(u,v,angle):
    # np.cos(angle*np.pi/180)
    # np.sin(angle*np.pi/180)
    u_east=u*np.cos(angle*np.pi/180)-v*np.sin(angle*np.pi/180)
    v_north=u*np.sin(angle*np.pi/180)+v*np.cos(angle*np.pi/180)
    return u_east,v_north


# In[ ]:





# In[ ]:





# In[18]:


grid_dir='/meddy/simingzhang/Data/RB_iceland_data/'
wave_dir='/meddy/simingzhang/Data/Parcels_data/'
# nowave_dir='/meddy/simingzhang/Data/RB_iceland_data/iceland_no_wave/'
grdname='niskin2km_500m_grd.nc'
# hisname_w='z_niskin2km_his_hf_depth_500m_grd.0002.nc'
# hisname_nw='z_niskin2km_his_smooth_depth_500m_grd.0002.nc'

grdname=f'{grid_dir}{grdname}'
# wavename=f'{wave_dir}{hisname_w}'
# nowavename=f'{nowave_dir}{hisname_nw}'

# wavename=f'{wave_dir}wavecase_modified_cg.nc'
# nowavename=f'{wave_dir}nowavecase_modified_cg.nc'
wavename=f'{wave_dir}wavecase_modified_vel_cg_tukey.nc'
nowavename=f'{wave_dir}nowavecase_modified_vel_cg_tukey.nc'


# In[ ]:





# In[19]:


grid=xr.open_dataset(grdname)
grid = grid.swap_dims({'eta_u': 'eta_rho','xi_v':'xi_rho'})

lon_rho=grid['lon_rho']
lat_rho=grid['lat_rho']
h=grid['h']
angle=grid['angle'].values
lonmin=np.min(lon_rho.values)
lonmax=np.max(lon_rho.values)
latmin=np.min(lat_rho.values)
latmax=np.max(lat_rho.values)
lonmin,lonmax,latmin,latmax
grid


# In[20]:


ds = xr.open_dataset(wavename)
ds
# depth=ds['depth'].values


# In[21]:


# u1=u2rho_3d(np.squeeze(ds['u'].values))
# v1=v2rho_3d(np.squeeze(ds['v'].values))
# ueast,vnorth=trans_vel_roms(u1,v1,angle)
# ds['u_rho'].values[:,0,:,:]=ueast
# ds['v_rho'].values[:,0,:,:]=vnorth
# ds.to_netcdf(f'{wave_dir}wavecase_modified_cg_uni1.nc')


# # Parcel
# ### 1. wave case

# In[22]:


wavename
# filenames = {"U": wavename, "V": wavename,"PI2": wavename,"PI4": wavename,"PI6": wavename,"PI8": wavename,"PI10": wavename,
#             "PI12": wavename,"PI16": wavename,"PI20": wavename,"PI30": wavename,"PI50": wavename,"PI60": wavename,"PI100": wavename}

# filenames = {"U": wavename, "V": wavename}

# filenames = {"U": wavename, "V": wavename,"Th1": wavename,"Th2": wavename,"Th3": wavename,"Th4": wavename,"Th5": wavename,
#             "Th6": wavename,"Th7": wavename,"Th8": wavename,"Th9": wavename,"Th10": wavename,"Th11": wavename,"Th12": wavename,
#             "Th13": wavename,"Th14": wavename,"Th15": wavename,"Th16": wavename,"Th17": wavename,"Th18": wavename,"Th21": wavename,
#             "Th24": wavename,"Th27": wavename,"Th30": wavename,"Th33": wavename,"Th36": wavename,"Th39": wavename,"Th42": wavename,
#              "Th45": wavename,"Th48": wavename,"Th54": wavename,"Th60": wavename,"Th66": wavename,"Th72": wavename,"Th78": wavename,
#              "Th84": wavename,"Th90": wavename,"Th96": wavename,"Th102": wavename,"Th108": wavename,"Th114": wavename}

# 创建变量名列表
variable_names = ["U", "V", "Ro", "divof"]
variable_names += [f"Th{i}" for i in range(1, 19)]
variable_names += [f"Th{i}" for i in [21, 24, 27, 30, 33, 36, 39, 42, 45, 48,
                                      54, 60, 66, 72, 78, 84, 90, 96, 102, 108, 114]]

# 创建字典
filenames = {var: wavename for var in variable_names}

# 打印结果
# print(filenames)
filenames


# In[ ]:





# In[ ]:





# In[23]:


# pip show parcels


# In[24]:


# variables = {
#     "U": "u_rho",
#     "V": "v_rho",
#     "PI2": "PI2",
#     "PI4": "PI4",
#     "PI6": "PI6",
#     "PI8": "PI8",
#     "PI10": "PI10",
#     "PI12": "PI12",
#     "PI16": "PI16",
#     "PI20": "PI20",
#     "PI30": "PI30",
#     "PI50": "PI50",
#     "PI60": "PI60",
#     "PI100": "PI100",
# }
# dimensions = {
#     "U": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "V": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "PI2": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "PI4": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "PI6": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "PI8": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "PI10": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "PI12": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "PI16": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "PI20": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "PI30": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "PI50": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "PI60": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "PI100": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
# }


# variables = {
#     "U": "u_rho","V": "v_rho",
#     "Th1": "Th1","Th2": "Th2","Th3": "Th3","Th4": "Th4","Th5": "Th5","Th6": "Th6","Th7": "Th7","Th8": "Th8","Th9": "Th9",
#     "Th10": "Th10","Th11": "Th11","Th12": "Th12","Th13": "Th13","Th14": "Th14","Th15": "Th15","Th16": "Th16","Th17": "Th17","Th18": "Th18",
#     "Th21": "Th21","Th24": "Th24","Th27": "Th27","Th30": "Th30","Th33": "Th33","Th36": "Th36","Th39": "Th39","Th42": "Th42","Th45": "Th45",
#     "Th48": "Th48","Th54": "Th54","Th60": "Th60","Th66": "Th66","Th72": "Th72","Th78": "Th78","Th84": "Th84","Th90": "Th90","Th96": "Th96",
#     "Th102": "Th102","Th108": "Th108","Th114": "Th114",
    
# }
# dimensions = {
#     "U": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "V": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th1": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th2": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th3": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th4": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th5": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th6": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th7": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th8": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th9": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th10": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th11": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th12": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th13": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th14": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th15": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th16": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th17": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th18": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th21": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th24": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th27": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th30": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th33": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th36": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th39": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th42": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th45": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th48": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th54": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th60": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th66": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th72": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th78": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th84": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th90": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th96": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th102": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th108": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},
#     "Th114": {"lat": "lat_rho", "lon": "lon_rho", "time": "ocean_time", "depth": "depth",},

# }


all_vars = ["U", "V", "Ro", "divof"] + [f"Th{i}" for i in list(range(1, 19)) + [21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 54, 60, 66, 72, 78, 84, 90, 96, 102, 108, 114]]

# 创建 variables 字典
variables = {
    var: "u_rho" if var == "U" else "v_rho" if var == "V" else "Ro" if var == "Ro" else "divof" if var == "divof" else var
    for var in all_vars
}

# 维度模板
dim_template = {
    "lat": "lat_rho",
    "lon": "lon_rho",
    "time": "ocean_time",
    "depth": "depth"
}

# 创建 dimensions 字典
dimensions = {var: dim_template.copy() for var in all_vars}

chunks=500
cs = {"time": ("ocean_time", 1), "lat": ("lat_rho", chunks), "lon": ("lon_rho", chunks),"depth": ("depth", chunks)}
# ,chunksize=cs
fieldset = parcels.FieldSet.from_netcdf(filenames, variables, dimensions,allow_time_extrapolation=True)
# ,allow_time_extrapolation=True
# 


# In[25]:


fieldset


# In[26]:


# fieldset.U.grid.time *= 3600  # 秒 → 小时
# lonp1.shape


# In[27]:


# fieldset.U.grid.time/1e9


# In[28]:


def DeleteParticle(particle, fieldset, time):
    if (particle.lon < -30.92595646869225+1.5) or (particle.lon > -19.9541352315797-1.5) or \
       (particle.lat < 56.70597152564239+0.5) or (particle.lat > 62.267609501148414-0.5):
        particle.delete()
        
# def SampleT(particle, fieldset, time):
#     particle.ue = fieldset.U[time, particle.depth, particle.lat, particle.lon]
#     particle.ve = fieldset.V[time, particle.depth, particle.lat, particle.lon]
#     particle.pi2 = fieldset.PI2[time, particle.depth, particle.lat, particle.lon]
#     particle.pi4 = fieldset.PI4[time, particle.depth, particle.lat, particle.lon]
#     particle.pi6 = fieldset.PI6[time, particle.depth, particle.lat, particle.lon]
#     particle.pi8 = fieldset.PI8[time, particle.depth, particle.lat, particle.lon]
#     particle.pi10 = fieldset.PI10[time, particle.depth, particle.lat, particle.lon]
#     particle.pi12 = fieldset.PI12[time, particle.depth, particle.lat, particle.lon]
#     particle.pi16 = fieldset.PI16[time, particle.depth, particle.lat, particle.lon]
#     particle.pi20 = fieldset.PI20[time, particle.depth, particle.lat, particle.lon]
#     particle.pi30 = fieldset.PI30[time, particle.depth, particle.lat, particle.lon]
#     particle.pi50 = fieldset.PI50[time, particle.depth, particle.lat, particle.lon]
#     particle.pi60 = fieldset.PI60[time, particle.depth, particle.lat, particle.lon]
#     particle.pi100 = fieldset.PI100[time, particle.depth, particle.lat, particle.lon]

def SampleT(particle, fieldset, time):
    particle.ue = fieldset.U[time, particle.depth, particle.lat, particle.lon]
    particle.ve = fieldset.V[time, particle.depth, particle.lat, particle.lon]
    particle.th1 = fieldset.Th1[time, particle.depth, particle.lat, particle.lon]
    particle.th2 = fieldset.Th2[time, particle.depth, particle.lat, particle.lon]
    particle.th3 = fieldset.Th3[time, particle.depth, particle.lat, particle.lon]
    particle.th4 = fieldset.Th4[time, particle.depth, particle.lat, particle.lon]
    particle.th5 = fieldset.Th5[time, particle.depth, particle.lat, particle.lon]
    particle.th6 = fieldset.Th6[time, particle.depth, particle.lat, particle.lon]
    particle.th7 = fieldset.Th7[time, particle.depth, particle.lat, particle.lon]
    particle.th8 = fieldset.Th8[time, particle.depth, particle.lat, particle.lon]
    particle.th9 = fieldset.Th9[time, particle.depth, particle.lat, particle.lon]
    particle.th10 = fieldset.Th10[time, particle.depth, particle.lat, particle.lon]
    particle.th11 = fieldset.Th11[time, particle.depth, particle.lat, particle.lon]
    particle.th12 = fieldset.Th12[time, particle.depth, particle.lat, particle.lon]
    particle.th13 = fieldset.Th13[time, particle.depth, particle.lat, particle.lon]
    particle.th14 = fieldset.Th14[time, particle.depth, particle.lat, particle.lon]
    particle.th15 = fieldset.Th15[time, particle.depth, particle.lat, particle.lon]
    particle.th16 = fieldset.Th16[time, particle.depth, particle.lat, particle.lon]
    particle.th17 = fieldset.Th17[time, particle.depth, particle.lat, particle.lon]
    particle.th18 = fieldset.Th18[time, particle.depth, particle.lat, particle.lon]
    particle.th21 = fieldset.Th21[time, particle.depth, particle.lat, particle.lon]
    particle.th24 = fieldset.Th24[time, particle.depth, particle.lat, particle.lon]
    particle.th27 = fieldset.Th27[time, particle.depth, particle.lat, particle.lon]
    particle.th30 = fieldset.Th30[time, particle.depth, particle.lat, particle.lon]
    particle.th33 = fieldset.Th33[time, particle.depth, particle.lat, particle.lon]
    particle.th36 = fieldset.Th36[time, particle.depth, particle.lat, particle.lon]
    particle.th39 = fieldset.Th39[time, particle.depth, particle.lat, particle.lon]
    particle.th42 = fieldset.Th42[time, particle.depth, particle.lat, particle.lon]
    particle.th45 = fieldset.Th45[time, particle.depth, particle.lat, particle.lon]
    particle.th48 = fieldset.Th48[time, particle.depth, particle.lat, particle.lon]
    particle.th54 = fieldset.Th54[time, particle.depth, particle.lat, particle.lon]
    particle.th60 = fieldset.Th60[time, particle.depth, particle.lat, particle.lon]
    particle.th66 = fieldset.Th66[time, particle.depth, particle.lat, particle.lon]
    particle.th72 = fieldset.Th72[time, particle.depth, particle.lat, particle.lon]
    particle.th78 = fieldset.Th78[time, particle.depth, particle.lat, particle.lon]
    particle.th84 = fieldset.Th84[time, particle.depth, particle.lat, particle.lon]
    particle.th90 = fieldset.Th90[time, particle.depth, particle.lat, particle.lon]
    particle.th96 = fieldset.Th96[time, particle.depth, particle.lat, particle.lon]
    particle.th102 = fieldset.Th102[time, particle.depth, particle.lat, particle.lon]
    particle.th108 = fieldset.Th108[time, particle.depth, particle.lat, particle.lon]
    particle.th114 = fieldset.Th114[time, particle.depth, particle.lat, particle.lon]
    particle.Ro = fieldset.Ro[time, particle.depth, particle.lat, particle.lon]
    particle.divof = fieldset.divof[time, particle.depth, particle.lat, particle.lon]



def CheckOutOfBounds(particle, fieldset, time):
    if particle.state == StatusCode.ErrorOutOfBounds:
        particle.delete()
def CheckError(particle, fieldset, time):
    if particle.state >= 50:  # This captures all Errors
        particle.delete()


# In[29]:


lonp=lon_rho.values
latp=lat_rho.values
Range=20
JJ=[Range,287-Range]
II=[Range,287-Range]
step=15
lonp1=lonp[JJ[0]:JJ[1]:step,II[0]:II[1]:step].ravel()
latp1=latp[JJ[0]:JJ[1]:step,II[0]:II[1]:step].ravel()
# step=1
# centerpoint=144
# # npart=2500
# npart=17*17
# nrelease=40
# repeatdt=6
# if (np.sqrt(npart) / 2) % 1 == 0:
#     left1=int(np.sqrt(npart)/2)
#     right1=int(np.sqrt(npart)/2)
#     # print(1)
# else:
#     left1=int((np.sqrt(npart)+1)/2)
#     right1=int((np.sqrt(npart)+1)/2)-1
#     # print(2)
# II=[centerpoint-left1,centerpoint+right1]
# JJ=[centerpoint-left1,centerpoint+right1]
lonp1=lonp[JJ[0]:JJ[1]:step,II[0]:II[1]:step].ravel()
latp1=latp[JJ[0]:JJ[1]:step,II[0]:II[1]:step].ravel()
during=89.5 # days
total_seconds = during * 86400 
# outputwavename=f'{wave_dir}wave_pars_P{latp1.shape[0]}T{during}days.zarr'
# lonp2= np.tile(lonp1, int(nrelease))
# latp2= np.tile(latp1, int(nrelease))
# timep2=np.repeat(np.arange(0, int(nrelease))*timedelta(hours=repeatdt).total_seconds(),npart)
# outputwavename=f'{wave_dir}wave_pars_P{int(lonp2.shape[0])}T{during}days.zarr'
outputwavename=f'{wave_dir}wave_pars_P{int(lonp1.shape[0])}T{during}days.zarr'


# In[30]:


lonp1.shape


# In[31]:


# lonp2.shape,timep2.shape


# In[32]:


# timep2


# In[33]:


# pset = parcels.ParticleSet.from_list(
#     fieldset=fieldset,
#     pclass=parcels.JITParticle.add_variable("ue").add_variable("ve").add_variable("th1").add_variable("th2").add_variable("th3").
#     add_variable("th4").add_variable("th5").add_variable("th6").add_variable("th7").add_variable("th8").
#     add_variable("th9").add_variable("th10").add_variable("th11").add_variable("th12").add_variable("th13").
#     add_variable("th14").add_variable("th15").add_variable("th16").add_variable("th17").add_variable("th18").
#     add_variable("th21").add_variable("th24").add_variable("th27").add_variable("th30").add_variable("th33").
#     add_variable("th36").add_variable("th39").add_variable("th42").add_variable("th45").add_variable("th48").
#     add_variable("th54").add_variable("th60").add_variable("th66").add_variable("th72").add_variable("th78").
#     add_variable("th84").add_variable("th90").add_variable("th96").add_variable("th102").add_variable("th108").
#     add_variable("th114"),
#     lon=lonp1,
#     lat=latp1,
#     # size=100,                # 粒子总数
#     depth=-2*np.ones(lonp1.shape)                  # 释放深度
# )

# output_file = pset.ParticleFile(
#     name=outputwavename, outputdt=timedelta(hours=1)
# )

# pset.execute(
#     [SampleT,parcels.AdvectionRK4,CheckOutOfBounds],  # simply combine the Kernels in a list
#     runtime=timedelta(hours=2148),
#     dt=timedelta(seconds=600),
#     output_file=output_file,
# )

# late start
assert len(lonp1) == len(latp1), "经度和纬度数组长度不一致"

# repeatdt = timedelta(hours=6)
pset = parcels.ParticleSet.from_list(
    fieldset=fieldset,
    pclass=parcels.JITParticle.add_variable("ue").add_variable("ve").add_variable("th1").add_variable("th2").add_variable("th3").
    add_variable("th4").add_variable("th5").add_variable("th6").add_variable("th7").add_variable("th8").
    add_variable("th9").add_variable("th10").add_variable("th11").add_variable("th12").add_variable("th13").
    add_variable("th14").add_variable("th15").add_variable("th16").add_variable("th17").add_variable("th18").
    add_variable("th21").add_variable("th24").add_variable("th27").add_variable("th30").add_variable("th33").
    add_variable("th36").add_variable("th39").add_variable("th42").add_variable("th45").add_variable("th48").
    add_variable("th54").add_variable("th60").add_variable("th66").add_variable("th72").add_variable("th78").
    add_variable("th84").add_variable("th90").add_variable("th96").add_variable("th102").add_variable("th108").
    add_variable("th114").add_variable("Ro").add_variable("divof"),
    lon=lonp1,
    lat=latp1,
    # size=100,                # 粒子总数
    depth=-2*np.ones(lonp1.shape),
    # time=timep2,
    # repeatdt=repeatdt,
    # chunks=(10000, 1),
)


output_file = pset.ParticleFile(
    name=outputwavename, outputdt=timedelta(hours=1)
)

pset.execute(
    [SampleT,parcels.AdvectionRK4,CheckOutOfBounds],  # simply combine the Kernels in a list
    runtime=timedelta(hours=2148),
    dt=timedelta(seconds=600),
    output_file=output_file,
)

# pset.repeatdt = None


# pset.execute(
#     [SampleT,parcels.AdvectionRK4,CheckOutOfBounds],  # simply combine the Kernels in a list
#     # runtime=timedelta(hours=9),
#     runtime=timedelta(hours=2098),
#     dt=timedelta(seconds=600),
#     output_file=output_file,
# )


# In[ ]:


# outputwavename=f'{wave_dir}wave_pars_P{100000}T{89.5}days.zarr'


# In[34]:


outputwavename


# In[35]:


# zarr_path = '/meddy/simingzhang/Data/Parcels_data/wave_pars_P103462.0T89.5days.zarr'

# # 使用xarray的open_mfdataset打开分布式存储
# ds = xr.open_mfdataset(
#     f'{zarr_path}/proc*.zarr',  # 匹配所有进程目录
#     engine='zarr',
#     combine='nested',
#     concat_dim='particle',
#     parallel=True
# )


# from glob import glob
# from os import path

# files = glob(path.join(f'{outputwavename}/', "proc*"))
# ds = xr.concat(
#     [xr.open_zarr(f) for f in files],
#     dim="trajectory",
#     compat="no_conflicts",
#     coords="minimal",
# )


# In[36]:


ds = xr.open_zarr(outputwavename)


# In[37]:


ds


# In[38]:


# plt.plot(ds.lon.T, ds.lat.T, ".-")
# plt.xlabel("lon")
# plt.ylabel("lat")
# plt.show()


# # transform velocity

# In[39]:


ue,ve=transunit_spher2flat(ds['ue'].values,ds['ve'].values,ds['lat'].values)
ds['ue'].values=ue
ds['ve'].values=ve


# In[40]:


f'{outputwavename[37:-5]}.nc'


# In[41]:


ds.to_netcdf(f'{outputwavename[:-5]}.nc')


# In[ ]:


# !jupyter nbconvert --to python Parcels_Iceland_wave.ipynb


# In[ ]:


# (ds['time'].values/1e9/3600)[:,0].astype(int)


# In[ ]:


# list(ds.data_vars)[10]


# In[ ]:


# ds['lon'].values.shape[0]


# In[ ]:


# LONt=np.zeros((ds['lon'].values.shape))
# for ii in range(ds['lon'].values.shape[0]):
#     TT=ds['time'].values[ii,:].astype(int)
#     indice=int(TT[0]/3600/1e9)
#     midvar=np.zeros(TT.shape)*np.nan;
#     midvar[0+indice:]=(ds['lon'].values)[ii,0:(TT.shape[0]-indice)];
#     LONt[ii,:]=midvar;
#     print(ii)
# ds['lon'].values=LONt


# In[26]:


# ds[list(ds.data_vars)[0]]


# In[27]:


# ds1=ds


# In[28]:


# ds1

