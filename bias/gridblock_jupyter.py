"""Python ports of Homo_gridblock_Iceland.m and homo_gridblock_hit.m."""
from pathlib import Path
import time
import numpy as np
import xarray as xr
from scipy.interpolate import LinearNDInterpolator
from scipy.io import savemat
from tqdm import tqdm

EARTH_RADIUS = 6367442.76


def matlab_layout(a, particle_count, name):
    a = np.asarray(a).squeeze()
    if a.ndim != 2:
        raise ValueError(f"{name} must be 2-D, got {a.shape}")
    if a.shape[1] == particle_count:
        return a
    if a.shape[0] == particle_count:
        return a.T
    raise ValueError(f"Cannot identify particle dimension for {name}: {a.shape}")


def spheric_dist(lat1, lat2, lon1, lon2):
    l = np.abs(lon2 - lon1)
    l = np.where(l >= 180, 360 - l, l)
    lat1, lat2, l = np.deg2rad(lat1), np.deg2rad(lat2), np.deg2rad(l)
    q = (np.sin(l)*np.cos(lat2))**2 + (np.sin(lat2)*np.cos(lat1) - np.sin(lat1)*np.cos(lat2)*np.cos(l))**2
    return EARTH_RADIUS*np.arcsin(np.sqrt(np.clip(q, 0, 1)))


def condensed_pairs(x, y, u, v, geographic):
    ii, jj = np.triu_indices(x.size, 1)
    rx = x[ii] - x[jj]
    if geographic:
        rx *= np.cos(np.deg2rad(0.5*(y[ii] + y[jj])))
        factor = 111321.0
    else:
        factor = 1.0
    ry = y[ii] - y[jj]
    magr = np.sqrt(rx*rx + ry*ry)
    dist = factor*magr
    dux, duy = u[ii] - u[jj], v[ii] - v[jj]
    with np.errstate(divide="ignore", invalid="ignore"):
        dul = dux*rx/magr + duy*ry/magr
        dut = duy*rx/magr - dux*ry/magr
    good = np.isfinite(dist) & np.isfinite(dul) & np.isfinite(dut) & np.isfinite(magr) & (magr > 0)
    return ii[good], jj[good], dist[good], dul[good], dut[good]


def discretize(values, edges):
    """Zero-based equivalent of MATLAB discretize; -1 means outside/NaN."""
    values = np.asarray(values)
    out = np.searchsorted(edges, values, side="right") - 1
    out[values == edges[-1]] = len(edges) - 2
    out[(values < edges[0]) | (values > edges[-1]) | ~np.isfinite(values)] = -1
    return out


def calculate_H(local_sf, pair_count, min_pairs, min_valid_blocks):
    ns = local_sf.shape[0]
    H = np.full(ns, np.nan); mean = np.full(ns, np.nan); std = np.full(ns, np.nan); rms = np.full(ns, np.nan)
    nvalid = np.zeros(ns, dtype=np.int64)
    for ir in range(ns):
        valid = (pair_count[ir] >= min_pairs) & np.isfinite(local_sf[ir])
        values = local_sf[ir, valid]
        nvalid[ir] = values.size
        if values.size < min_valid_blocks:
            continue
        mean[ir] = values.mean()
        std[ir] = values.std(ddof=0)  # MATLAB std(values,1)
        rms[ir] = np.sqrt(np.mean(values**2))
        denominator = abs(mean[ir]) + rms[ir]
        if np.isfinite(denominator) and denominator > 0:
            H[ir] = std[ir]/denominator
    return H, mean, std, rms, nvalid


def finish_statistics(sums, count_total, block_sums, pair_count, min_pairs, min_valid_blocks):
    overall = []
    valid = count_total > 0
    for total in sums:
        a = np.full(count_total.size, np.nan); a[valid] = total[valid]/count_total[valid]; overall.append(a)
    local = []
    valid_local = pair_count >= min_pairs
    for total in block_sums:
        a = np.full(pair_count.shape, np.nan); a[valid_local] = total[valid_local]/pair_count[valid_local]; local.append(a)
    hstats = [calculate_H(a, pair_count, min_pairs, min_valid_blocks) for a in local]
    return overall, local, hstats


def add_by_bin(target, bins, values):
    good = bins >= 0
    np.add.at(target, bins[good], values[good])


def matlab_bins(first_edge, gamma, cutoff, zero_first=True):
    raw = first_edge*gamma**np.arange(101)
    first_large = np.flatnonzero(raw > cutoff)[0]
    edges = np.r_[0.0, raw[:first_large]] if zero_first else np.r_[raw[2]-raw[1], raw[:first_large]]
    return edges, 0.5*(edges[:-1] + edges[1:])


def bootstrap_block_mean(samples, num_boot, rng, chunk=64):
    nblock = samples.shape[0]
    out = np.empty(num_boot)
    for start in range(0, num_boot, chunk):
        stop = min(start + chunk, num_boot)
        draw = rng.integers(0, nblock, size=(stop-start, nblock))
        out[start:stop] = samples[draw].mean(axis=(1, 2))
    return out


def run_iceland(cfg):
    t0 = time.perf_counter(); rng = np.random.default_rng(cfg.get("random_seed"))
    tr = np.asarray(cfg["timerange_matlab"], int)-1; npart = cfg["nparticles"]
    fname = Path(cfg["input_dir"])/f'{cfg["case"]}_pars_P{npart}T{cfg["days"]:g}days.nc'
    with xr.open_dataset(fname, decode_times=False) as ds:
        lon0 = matlab_layout(ds.lon.values, npart, "lon")[tr]; lat0 = matlab_layout(ds.lat.values, npart, "lat")[tr]
        _ = matlab_layout(ds.ue.values, npart, "ue")[tr]; _ = matlab_layout(ds.ve.values, npart, "ve")[tr]
    dt = cfg["dt"]
    U = spheric_dist(lat0[:-1], lat0[:-1], lon0[:-1], lon0[1:])/dt*np.sign(lon0[1:]-lon0[:-1])
    V = spheric_dist(lat0[1:], lat0[:-1], lon0[:-1], lon0[:-1])/dt*np.sign(lat0[1:]-lat0[:-1])
    lon, lat = lon0[:-1], lat0[:-1]; ntime = lon.shape[0]
    with xr.open_dataset(cfg["grid_file"], decode_times=False) as grid:
        # MATLAB ncread reverses NetCDF dimension order; transpose to reproduce it.
        lr = np.asarray(grid.lon_rho.values).T; ar = np.asarray(grid.lat_rho.values).T
    ni, nj = lr.shape; ig, jg = np.meshgrid(np.arange(1, ni+1), np.arange(1, nj+1), indexing="ij")
    vg = np.isfinite(lr) & np.isfinite(ar); pts = np.c_[lr[vg], ar[vg]]
    FI = LinearNDInterpolator(pts, ig[vg], fill_value=np.nan); FJ = LinearNDInterpolator(pts, jg[vg], fill_value=np.nan)
    vp = np.isfinite(lon) & np.isfinite(lat); pp = np.c_[lon[vp], lat[vp]]; ia, ja = FI(pp), FJ(pp); ok = np.isfinite(ia) & np.isfinite(ja)
    if not ok.any(): raise ValueError("Particle positions cannot be mapped to model grid")
    imin, imax = max(1, np.floor(ia[ok].min())), min(ni, np.ceil(ia[ok].max()))
    jmin, jmax = max(1, np.floor(ja[ok].min())), min(nj, np.ceil(ja[ok].max()))
    ie = np.linspace(imin, np.nextafter(imax, np.inf), cfg["nblock_I"]+1); je = np.linspace(jmin, np.nextafter(jmax, np.inf), cfg["nblock_J"]+1)
    r_requested_km = np.asarray(cfg["r_requested_km"], float).ravel()
    r_requested_km = np.sort(np.unique(r_requested_km[np.isfinite(r_requested_km) & (r_requested_km > 0)]))
    if r_requested_km.size < 2:
        raise ValueError("r_requested_km must contain at least two positive scales")
    lr = np.log(r_requested_km)
    le = np.empty(lr.size + 1); le[1:-1] = 0.5*(lr[:-1] + lr[1:])
    le[0] = lr[0] - 0.5*(lr[1]-lr[0]); le[-1] = lr[-1] + 0.5*(lr[-1]-lr[-2])
    db, da = np.exp(le)*1000.0, r_requested_km*1000.0
    ns, nb = da.size, cfg["nblock_I"]*cfg["nblock_J"]
    sums = [np.zeros(ns) for _ in range(5)]; bsums = [np.zeros((ns, nb)) for _ in range(3)]; counts = np.zeros(ns); pcount = np.zeros((ns, nb))
    pairs = [[[], []] for _ in range(ns)] if cfg["do_bootstrap"] else None
    for it in tqdm(
        range(ntime),
        total=ntime,
        desc="Iceland Lagrangian pair statistics",
        unit="time"
    ):
        valid = np.isfinite(lon[it]) & np.isfinite(lat[it]) & np.isfinite(U[it]) & np.isfinite(V[it]); ids = np.flatnonzero(valid)
        if ids.size < 2: continue
        x,y,u,v = lon[it,ids],lat[it,ids],U[it,ids],V[it,ids]; ii,jj,d,dl,dtr = condensed_pairs(x,y,u,v,True)
        mid = np.c_[.5*(x[ii]+x[jj]), .5*(y[ii]+y[jj])]; si = discretize(d,db); bi = discretize(FI(mid),ie); bj = discretize(FJ(mid),je)
        block = bi*cfg["nblock_J"]+bj; good=(si>=0)&(bi>=0)&(bj>=0); si,block,dl,dtr=si[good],block[good],dl[good],dtr[good]
        vals=[dl,dl**2,dtr**2,dl**3,dl*dtr**2]
        for a,z in zip(sums,vals): add_by_bin(a,si,z)
        add_by_bin(counts,si,np.ones(si.size))
        flat=si*nb+block
        for a,z in zip(bsums,[dl,dl**2,dl**3]): np.add.at(a.ravel(),flat,z)
        np.add.at(pcount.ravel(),flat,1)
        if pairs is not None:
            for ir in np.unique(si):
                q=si==ir; pairs[ir][0].append(dl[q]); pairs[ir][1].append(dtr[q])
    pair_runtime=time.perf_counter()-t0
    overall, local, hs = finish_statistics(sums,counts,bsums,pcount,cfg["min_pairs"],cfg["min_valid_blocks"])
    SF1,SF2ll,SF2tt,SF3lll,SF3ltt=overall; SF2=SF2ll+SF2tt; SF3=SF3lll+SF3ltt
    out={"Case":cfg["case"],"nparticles":npart,"days":cfg["days"],"dt":dt,"timerange":tr+1,"r_requested_km":r_requested_km,"dist_axis":da/1000.0,"dist_bin":db/1000.0,
         "SF1":SF1,"SF2":SF2,"SF2ll":SF2ll,"SF2tt":SF2tt,"SF3":SF3,"SF3lll":SF3lll,"SF3ltt":SF3ltt,"count_total":counts,
         "local_SF1":local[0],"local_SF2":local[1],"local_SF3":local[2],"pair_count":pcount,"I_edges":ie,"J_edges":je,"I_min":imin,"I_max":imax,"J_min":jmin,"J_max":jmax,
         "nblock_I":cfg["nblock_I"],"nblock_J":cfg["nblock_J"],"nBlock":nb,"min_pairs":cfg["min_pairs"],"min_valid_blocks":cfg["min_valid_blocks"],"do_bootstrap":cfg["do_bootstrap"],"num_boot":cfg["num_boot"],"pair_runtime":pair_runtime}
    for k,h in zip(("1","2","3"),hs):
        for name,val in zip(("H","mean_SF","std_SF","rms_SF","nvalid"),h): out[name+k]=val
    if pairs is not None:
        bt=time.perf_counter(); nboot=cfg["num_boot"]; matrices=[np.full((ns,nboot),np.nan) for _ in range(4)]; dof=np.ones(ns); good=np.isfinite(SF2)&(SF2>0); dof[good]=np.ceil((ntime*dt)/(da[good]/np.sqrt(SF2[good]))); dof=np.maximum(dof,1); sample=np.zeros(ns)
        for ir,p in enumerate(pairs):
            dl=np.concatenate(p[0]) if p[0] else np.array([]); dtr=np.concatenate(p[1]) if p[1] else np.array([]); sample[ir]=dl.size
            if dl.size<=10: continue
            nblocks=int(min(dof[ir],dl.size)); blocksize=dl.size//nblocks; order=rng.choice(dl.size,nblocks*blocksize,replace=False); bd=dl[order].reshape(nblocks,blocksize); bt_=dtr[order].reshape(nblocks,blocksize)
            for mat,z in zip(matrices,(bd,bd**2,bd**3,bd**3+bd*bt_**2)): mat[ir]=bootstrap_block_mean(z,nboot,rng)
        out.update({"dof":dof,"SF1l":matrices[0],"SF2l":matrices[1],"SF3l":matrices[2],"SF3full_boot":matrices[3],"SF3_mean":np.nanmean(matrices[3],axis=1),"SF3_stderr":np.nanstd(matrices[3],axis=1,ddof=1),"nsample":sample,"bootstrap_runtime":time.perf_counter()-bt})
    else:
        out.update({k:np.array([]) for k in ("dof","SF1l","SF2l","SF3l","SF3full_boot","SF3_mean","SF3_stderr","nsample")}); out["bootstrap_runtime"]=np.nan
    path=Path(cfg["input_dir"])/f'{cfg["case"]}_pars_P{npart}T{tr[-1]+1}{cfg["ini"]}gridBlockHL.mat'; savemat(path,out,do_compression=True,oned_as="row"); print("saved",path); return out,path


def run_hit(cfg):
    t0=time.perf_counter(); rng=np.random.default_rng(cfg["random_seed"]); tr=np.asarray(cfg["timerange_matlab"],int)-1
    fname=Path(cfg["input_dir"])/f'HIT2d_pars_P{cfg["nparticles_file"]}T{cfg["seconds"]:.1f}seconds.nc'
    with xr.open_dataset(fname,decode_times=False) as ds:
        lon0=matlab_layout(ds.lon.values,cfg["nparticles_file"],"lon"); available=lon0.shape[1]
        selected=rng.choice(available,cfg["num_to_select"],replace=False)
        lon=lon0[np.ix_(tr,selected)]; lat=matlab_layout(ds.lat.values,available,"lat")[np.ix_(tr,selected)]
        u=matlab_layout(ds.ue.values,available,"ue")[np.ix_(tr,selected)]; v=matlab_layout(ds.ve.values,available,"ve")[np.ix_(tr,selected)]
        xs=np.asarray(cfg["xscale"],int); Th=np.vstack([np.nanmean(matlab_layout(ds[f"th{s}"].values,available,f"th{s}")[np.ix_(tr,selected)],axis=1) for s in xs])
    xd=np.array([np.nanmin(lon),np.nanmax(lon)]); yd=np.array([np.nanmin(lat),np.nanmax(lat)]); xe=np.linspace(*xd,cfg["nblock_x"]+1); ye=np.linspace(*yd,cfg["nblock_y"]+1)
    db,da=matlab_bins(.0123,1.2,3.14,False); ns,nb=da.size,cfg["nblock_x"]*cfg["nblock_y"]
    sums=[np.zeros(ns) for _ in range(5)]; bsums=[np.zeros((ns,nb)) for _ in range(3)]; counts=np.zeros(ns); pc=np.zeros((ns,nb)); particle_count=np.zeros((tr.size,nb))
    for it in range(tr.size):
        valid=np.isfinite(lon[it])&np.isfinite(lat[it])&np.isfinite(u[it])&np.isfinite(v[it]); bx=discretize(lon[it],xe); by=discretize(lat[it],ye); pb=bx*cfg["nblock_y"]+by
        for ib in range(nb):
            ids=np.flatnonzero(valid&(pb==ib));
            if ids.size<2: continue
            if ids.size>cfg["max_particles_per_block"]: ids=rng.choice(ids,cfg["max_particles_per_block"],replace=False)
            particle_count[it,ib]=ids.size; _,_,d,dl,dtr=condensed_pairs(lon[it,ids],lat[it,ids],u[it,ids],v[it,ids],False); si=discretize(d,db); good=si>=0; si,dl,dtr=si[good],dl[good],dtr[good]
            vals=[dl,dl**2,dtr**2,dl**3,dl*dtr**2]
            for a,z in zip(sums,vals): add_by_bin(a,si,z)
            add_by_bin(counts,si,np.ones(si.size)); flat=si*nb+ib
            for a,z in zip(bsums,[dl,dl**2,dl**3]): np.add.at(a.ravel(),flat,z)
            np.add.at(pc.ravel(),flat,1)
        print(f"time {it+1}/{tr.size}")
    runtime=time.perf_counter()-t0; overall,local,hs=finish_statistics(sums,counts,bsums,pc,cfg["min_pairs"],cfg["min_valid_blocks"]); SF1,SF2ll,SF2tt,SF3lll,SF3ltt=overall
    out={"Case":"HIT2d","nparticles_file":cfg["nparticles_file"],"num_to_select":cfg["num_to_select"],"selected_indices":selected+1,"seconds":cfg["seconds"],"dt":cfg["dt"],"timerange":tr+1,"xscale":xs,"Th_all":Th,"x_domain":xd,"y_domain":yd,"x_edges":xe,"y_edges":ye,"nblock_x":cfg["nblock_x"],"nblock_y":cfg["nblock_y"],"nBlock":nb,"max_particles_per_block":cfg["max_particles_per_block"],"min_pairs":cfg["min_pairs"],"min_valid_blocks":cfg["min_valid_blocks"],"dist_axis":da,"dist_bin":db,"SF1":SF1,"SF2":SF2ll+SF2tt,"SF2ll":SF2ll,"SF2tt":SF2tt,"SF3":SF3lll+SF3ltt,"SF3lll":SF3lll,"SF3ltt":SF3ltt,"local_SF1":local[0],"local_SF2":local[1],"local_SF3":local[2],"pair_count":pc,"count_total":counts,"particle_count":particle_count,"pair_runtime":runtime}
    for k,h in zip(("1","2","3"),hs):
        for name,val in zip(("H","mean_SF","std_SF","rms_SF","nvalid"),h): out[name+k]=val
    path=Path(cfg["input_dir"])/f'HIT2d_pars_P{cfg["num_to_select"]}T{tr[-1]+1}_gridBlockHL.mat'; savemat(path,out,do_compression=True,oned_as="row"); print("saved",path); return out,path


def max_curvature(log_res, log_sol):
    """Port of max_curvature.m using MATLAB-compatible gradient spacing."""
    dx = np.gradient(np.asarray(log_res, float)); dy = np.gradient(np.asarray(log_sol, float))
    d2x, d2y = np.gradient(dx), np.gradient(dy)
    with np.errstate(divide="ignore", invalid="ignore"):
        curvature = np.abs(d2y*dx-d2x*dy)/(dx*dx+dy*dy)**1.5
    return curvature, int(np.nanargmax(curvature))


def fit_sf3_lcurve(sf3, dist_axis, mindist, maxdist, lambda_vec, plot=False):
    """Port of the 'log' + 'RLS' path used by Fk_fitting_SF3_Lcurve.m."""
    from scipy.special import j1
    sf3=np.asarray(sf3,float); sf3=sf3[:,None] if sf3.ndim==1 else sf3
    r=np.asarray(dist_axis,float).ravel(); use=(r>=mindist)&(r<=maxdist); R=r[use]
    if R.size==0: raise ValueError("No valid points in requested distance range")
    kedge=np.logspace(np.log10(1/R.max()),np.log10(1/R.min()),R.size-1)*2*np.pi
    dk=np.diff(kedge)[::-1]; kf=(.5*(kedge[:-1]+kedge[1:]))[::-1]; nk=kf.size
    A=np.empty((R.size,nk+1))
    for j in range(nk): A[:,j]=-4/kf[j]*j1(kf[j]*R)*dk[j]
    A[:,-1]=2*R; A=A/np.abs(R)[:,None]
    ebs=np.zeros((nk+1,sf3.shape[1])); Vt=np.zeros((R.size,sf3.shape[1])); flux=np.zeros((nk,sf3.shape[1]))
    residuals=[]; solutions=[]; lam=np.asarray(lambda_vec,float)
    for value in lam:
        x=np.linalg.lstsq(np.vstack((A,np.sqrt(value)*np.eye(nk+1))),np.r_[sf3[use,0]/np.abs(R),np.zeros(nk+1)],rcond=None)[0]
        residuals.append(np.linalg.norm(A@x-sf3[use,0]/np.abs(R))); solutions.append(np.linalg.norm(x))
    _,idx=max_curvature(np.log10(residuals),np.log10(solutions)); optimal=lam[idx]
    for n in range(sf3.shape[1]):
        V=sf3[use,n]/np.abs(R); aug=np.vstack((A,np.sqrt(optimal)*np.eye(nk+1)))
        ebs[:,n]=np.linalg.lstsq(aug,np.r_[V,np.zeros(nk+1)],rcond=None)[0]
        Vt[:,n]=2*ebs[-1,n]*R
        for j in range(nk): Vt[:,n]-=4*ebs[j,n]/kf[j]*j1(kf[j]*R)*dk[j]
        flux[0,n]=-ebs[-1,n]
        for j in range(1,nk): flux[j,n]=flux[j-1,n]+ebs[-j,n]*dk[-j]
    flux=np.flipud(flux)
    if plot:
        import matplotlib.pyplot as plt
        plt.loglog(residuals,solutions,"-o"); plt.plot(residuals[idx],solutions[idx],"ro")
        for a,b,l in zip(residuals,solutions,lam): plt.text(a,b,f"lambda={l:.1e}")
        plt.xlabel("Residual norm"); plt.ylabel("Solution norm"); plt.grid(True)
    return flux,Vt,ebs,kf,R,optimal


def run_hit_timeavg2(cfg):
    """Port of the processing/join sections of process_parcels_HIT2d_timeavg2.m."""
    t0=time.perf_counter(); rng=np.random.default_rng(cfg.get("random_seed")); tr=np.asarray(cfg["timerange_matlab"],int)-1
    root=Path(cfg["input_dir"]); npart=cfg["nparticles"]; nselect=cfg["num_to_select"]
    fname=root/f'HIT2d_pars_P{npart}T{cfg["seconds"]:.1f}seconds.nc'
    with xr.open_dataset(fname,decode_times=False) as ds:
        lon0=matlab_layout(ds.lon.values,npart,"lon"); available=lon0.shape[1]
        selected=rng.choice(available,nselect,replace=False)
        ix=np.ix_(tr,selected); lons=lon0[ix]; lats=matlab_layout(ds.lat.values,available,"lat")[ix]
        ues=matlab_layout(ds.ue.values,available,"ue")[ix]; ves=matlab_layout(ds.ve.values,available,"ve")[ix]
        xs=np.asarray(cfg["xscale"],int); Th=np.vstack([np.nanmean(matlab_layout(ds[f"th{s}"].values,available,f"th{s}")[ix],axis=1) for s in xs])
    oname=root/f'HIT2d_pars_P{nselect}T{cfg["seconds"]:g}seconds.nc'
    savemat(oname.with_name(oname.stem+"traj.mat"),{"lons":lons,"lats":lats,"ues":ues,"ves":ves,"Th_all":Th},do_compression=True)
    H=np.full(ues.shape,-520.0); T_axis=np.linspace(cfg["dt"],ues.shape[0]*cfg["dt"],ues.shape[0])/86400
    traj={"trajmat_X":lons,"trajmat_Y":lats,"trajmat_U":ues,"trajmat_V":ves,"H":H,"T_axis":T_axis}
    savemat(fname.with_name(fname.stem+"traj.mat"),{"traj":traj},do_compression=True,oned_as="row")
    db,da=matlab_bins(.0123,1.2,3.14,False); ns=da.size; nt=tr.size
    fields=["SF1l","SF1t","SF2ll","SF2tt","SF3ltt","SF3lll"]
    time_values={k:np.full((ns,nt),np.nan) for k in fields}; nvalid_time=np.zeros((ns,nt))
    for it in range(nt):
        valid=H[it]<-500; ids=np.flatnonzero(valid)
        _,_,d,dl,dtr=condensed_pairs(lons[it,ids],lats[it,ids],ues[it,ids],ves[it,ids],False)
        bi=discretize(d,db); sums={k:np.zeros(ns) for k in fields}; nv=np.zeros(ns)
        values=(dl,dtr,dl**2,dtr**2,dl*dtr**2,dl**3)
        for key,z in zip(fields,values): add_by_bin(sums[key],bi,np.nan_to_num(z,nan=0.0))
        add_by_bin(nv,bi,np.ones(bi.size)); nvalid_time[:,it]=nv
        for key in fields:
            with np.errstate(divide="ignore",invalid="ignore"): time_values[key][:,it]=sums[key]/nv
        if cfg.get("save_chunks",False):
            c=it+1; cpath=fname.with_name(fname.stem+f"chunk{c:02d}traj.mat"); one={k:v[it:it+1] for k,v in traj.items() if k!="T_axis"}; one["T_axis"]=T_axis[it:it+1]
            savemat(cpath,{"traj":one},do_compression=True,oned_as="row")
            spath=fname.with_name(fname.stem+f"chunk{c:02d}SF123.mat"); savemat(spath,{**sums,"nvaild":nv,"dist_axis":da,"time_batch":1,"total_steps":nt},do_compression=True,oned_as="column")
        print(f"time {it+1}/{nt}")
    for a in time_values.values(): a[a==0]=np.nan
    SF3_mean=np.nanmean(time_values["SF3lll"]+time_values["SF3ltt"],axis=1)
    base={"SF3_mean":SF3_mean,"dist_axis":da,"Th_all":Th,"nvaild_time":nvalid_time}
    path5=root/f'HIT2d_pars_P{npart}T{tr[-1]+1}timeavg5.mat'; savemat(path5,base,do_compression=True,oned_as="row")
    path6=root/f'HIT2d_pars_P{npart}T{tr[-1]+1}timeavg6.mat'; savemat(path6,{**base,"SF1l_time":time_values["SF1l"],"SF2ll_time":time_values["SF2ll"],"SF3lll_time":time_values["SF3lll"]},do_compression=True,oned_as="row")
    print(f"finished in {time.perf_counter()-t0:.1f}s; saved {path5} and {path6}")
    return {**base,**{k+"_time":v for k,v in time_values.items()},"selected_indices":selected+1},(path5,path6)
