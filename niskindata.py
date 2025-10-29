
def path(grd, season, opt = 'hf'):
    if opt == 'hf':
        return hf(grd, season)
    elif opt=='smooth':
        return smooth(grd, season)
    elif opt=='hit':
        return hit2d(grd, season)
    elif opt=='aviso':
        return aviso(grd, season)
    elif opt=='hfrot':
        return hfrot(grd, season)
    elif opt=='smoothrot':
        return smoothrot(grd, season)

def smooth(grd, season):
    dr0 = '/meddy/simingzhang/Data/RB_iceland_data/'
#    dro0 = '/home/kaushiks/scratch/niskine/smooth/'
    dro0 = dr0
    if grd == '500':
            dr = '/meddy/simingzhang/Data/iceland_500m/'
            dr0 = '/meddy/simingzhang/Data/iceland_500m/'
            dro0 = dr0
            dro = dro0 + 'iceland_nowave/'        
            gridfile='sample_niskin_500m_grd.nc'
            filehis = 'iceland_nowave/afterspinup_z_his_depth'
            filez = 'z_uzvzbz_depth'
    else:
            dr = dr0         
            dro = dro0 + 'iceland_no_wave/'        
            gridfile='niskin2km_500m_grd.nc'
            if season == 'winter':
                filehis = 'z_niskin2km_his_smooth_depth_500m_grd'
                filez = 'z_uzvzbz_depth_500m_grd'
            else:
                filehis = 'z_niskin2km_his_smooth_summer_depth_500m_grd' 
                filez = 'z_uzvzbz_summer_depth_500m_grd'
    return dr, dro, gridfile, filehis, filez 

def hf(grd, season):
    dr0 = '/meddy/simingzhang/Data/RB_iceland_data/'
    dro0 = dr0
    #dro0 = '/home/kaushiks/scratch/niskin_flux/hf/'
    if grd == '500':
            dr = '/meddy/simingzhang/Data/iceland_500m/'
            dr0 = '/meddy/simingzhang/Data/iceland_500m/'
            dro0 = dr0
            dro = dro0 + 'iceland_wave/'        
            gridfile='sample_niskin_500m_grd.nc'
            if season=='summer':
                filehis = 'iceland_wave/z_niskin_500m_his_summer_depth_interp' 
                filez  = 'iceland_wave/z_uzvzbz_summer_depth_interp' 
            else:
                filehis = 'iceland_wave/afterspinup_z_his_depth'
                filez  = 'iceland_wave/z_uzvzbz_depth'
    elif grd == '2km':
            dr = dr0 
            dro = dro0 + 'iceland_wave/'       
            gridfile='niskin2km_500m_grd.nc'
            if season == 'winter':
                filehis = 'z_niskin2km_his_hf_depth_500m_grd'
                filez  = 'z_uzvzbz_depth_500m_grd'
            else:
                filehis = 'z_niskin2km_his_hf_summer_depth_interp_500m_grd'
                filez  = 'z_uzvzbz_summer_depth_interp_500m_grd'
    return dr, dro, gridfile, filehis, filez 

def hit2d(grd, season):
    dr0 = '/meddy/simingzhang/Data/RB_iceland_data/'
    dro0 = dr0
    #dro0 = '/home/kaushiks/scratch/niskin_flux/hf/'
    if grd == '500':
            dr = dr0 
            dro = dro0 + 'iceland_wave/'        
            gridfile='HIT2d_new.0002.nc'
            if season=='summer':
                filehis = 'HIT2d_new' 
                filez  = 'HIT2d_new' 
            else:
                filehis = 'HIT2d_new'
                filez  = 'HIT2d_new'
    elif grd == '2km':
            dr = dr0 
            dro = dro0 + 'iceland_wave/'       
            gridfile='HIT2d_new.0002.nc'
            if season == 'winter':
                filehis = 'HIT2d_new'
                filez  = 'HIT2d_new'
            else:
                filehis = 'HIT2d_new'
                filez  = 'HIT2d_new'
    return dr, dro, gridfile, filehis, filez 

def aviso(grd, season):
    dr0 = '/meddy/simingzhang/Data/RB_iceland_data/'
    dro0 = dr0
    #dro0 = '/home/kaushiks/scratch/niskin_flux/hf/'
    if grd == '500':
            dr = dr0 
            dro = dro0 + 'iceland_wave/'        
            gridfile='sample_niskin_500m_grd.nc'
            if season=='summer':
                filehis = 'AVISO_niskin2km' 
                filez  = 'AVISO_niskin2km' 
            else:
                filehis = 'AVISO_niskin2km'
                filez  = 'AVISO_niskin2km'
    elif grd == '2km':
            dr = dr0 
            dro = dro0 + 'iceland_wave/'       
            # gridfile='AVISO_niskin2km.0002.nc'
            gridfile='niskin2km_500m_grd.nc'
       
            if season == 'winter':
                filehis = 'AVISO_niskin2km'
                filez  = 'AVISO_niskin2km'
            else:
                filehis = 'AVISO_niskin2km'
                filez  = 'AVISO_niskin2km'
    return dr, dro, gridfile, filehis, filez 


def hfrot(grd, season):
    dr0 = '/meddy/simingzhang/Data/RB_iceland_data/'
#    dro0 = '/home/kaushiks/scratch/niskine/smooth/'
    dro0 = dr0
    if grd == '500':
            dr = '/meddy/simingzhang/Data/iceland_500m/'
            dr0 = '/meddy/simingzhang/Data/iceland_500m/'
            dro0 = dr0
            dro = dro0 + 'iceland_nowave/'        
            gridfile='sample_niskin_500m_grd.nc'
            filehis = 'iceland_nowave/afterspinup_z_his_depth'
            filez = 'z_uzvzbz_depth'
    else:
            dr = dr0         
            dro = dro0 + 'iceland_wave/'        
            gridfile='niskin2km_500m_grd.nc'
            if season == 'winter':
                filehis = 'iceland_wave/rot_helmholtz'
                filez = 'z_uzvzbz_depth_500m_grd'
            else:
                filehis = 'z_niskin2km_his_smooth_summer_depth_500m_grd' 
                filez = 'z_uzvzbz_summer_depth_500m_grd'
    return dr, dro, gridfile, filehis, filez 

def smoothrot(grd, season):
    dr0 = '/meddy/simingzhang/Data/RB_iceland_data/'
#    dro0 = '/home/kaushiks/scratch/niskine/smooth/'
    dro0 = dr0
    if grd == '500':
            dr = '/meddy/simingzhang/Data/iceland_500m/'
            dr0 = '/meddy/simingzhang/Data/iceland_500m/'
            dro0 = dr0
            dro = dro0 + 'iceland_nowave/'        
            gridfile='sample_niskin_500m_grd.nc'
            filehis = 'iceland_nowave/afterspinup_z_his_depth'
            filez = 'z_uzvzbz_depth'
    else:
            dr = dr0         
            dro = dro0 + 'iceland_no_wave/'        
            gridfile='niskin2km_500m_grd.nc'
            if season == 'winter':
                filehis = 'iceland_no_wave/rot_helmholtz'
                filez = 'z_uzvzbz_depth_500m_grd'
            else:
                filehis = 'z_niskin2km_his_smooth_summer_depth_500m_grd' 
                filez = 'z_uzvzbz_summer_depth_500m_grd'
    return dr, dro, gridfile, filehis, filez 
