% nparticles=[289,625,2500,15376]; % numbers of particles
days=89.5;  % days
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
input_dir='/meddy/simingzhang/Data/Parcels_data/';
inv_style='RLS';
rescale=1;
fitstr=1;
lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
    1e-7,1e-8,1e-9,1e-10,1e-11,1e-12];
N=287;dstr=1;dot=201;
colors={'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE'};
colors_rgb = {...
    [0, 114/255, 189/255], ...   % #0072BD
    [217/255, 83/255, 25/255], ... % #D95319
    [237/255, 177/255, 32/255], ... % #EDB120
    [126/255, 47/255, 142/255], ... % #7E2F8E
    [0.4660, 0.6740, 0.1880],...
    [0.3010, 0.7450, 0.9330]
};

if strcmpi(Case, 'wave')
    % cgname='s2sflux_spec_hf.0002.nc';
    cgname=['s2sflux_spec_hfrot.0002.nc'];
    % cgname=['s2sflux_spec_hf_',win,'.0002.nc'];
    % cgstd='wavecase_modified_cg_tukey1_RodivofstrainofCG_Eul_std_.mat'
end

if strcmpi(Case, 'nowave')
    cgname=['s2sflux_spec_smoothrot.0002.nc'];

    
    % cgname=['s2sflux_spec_smooth_',win,'.0002.nc'];
    % cgstd='nowavecase_modified_cg_tukey1_RodivofstrainofCG_Eul_std_.mat'
end

% cgname='s2sflux_spec_smooth.0002.nc';
ncdisp(cgname)
Thm_Eulerian=ncread(cgname,'Thm');
filtscale=ncread(cgname,'filtscale');
Thm_Eulerian=Thm_Eulerian(2:end);
filtscale=filtscale(2:end);



kf1=fliplr(1./filtscale'.*2.*pi);



eval(['load ',Case,'_rot_Eulerian_SF3.mat'])
% [r,SF3,S3L1,S3T1]=calc_radial(S3L1,S3T1,N,xscale);
S3L1=vertcat(S3L1_alltime{:});
S3T1=vertcat(S3T1_alltime{:});
r=r(2:end);
SF3=(S3T1(:,2:end)+S3L1(:,2:end))';
SF3_mean=nanmean(SF3,2);

% SF3=nanmean(S3L1_alltime(timerange,:)+S3T1_alltime(timerange,:),1);
clear SpecFlux_E;clear kf_E;
[SpecFlux_E,Vt_E,ebs_E,kf_E,lf_E]=Fk_fitting_SF3_Lcurve(SF3_mean(fitstr:dot), ...
    r(fitstr:dot),1,500e3,'fuc','RLS',lambda,kf1,0);

% [residual_norms, solution_norms,...
% lambda_opt_idx]=Fk_fitting_SF3_Lcurve(SF3(fitstr:dot)', ...
%     r(fitstr:dot),1,500e3, ...
% 'fuc','RLS',[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
% 1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14,1e-15],kf1);


kc=flipud(1./filtscale);
Thmc=flipud(Thm_Eulerian);
kc_mid_s=0.5.*(kc(2:end)+kc(1:end-1));
dPidk_s=(Thmc(2:end)-Thmc(1:end-1))

dk_E=u2rho_2d(abs(diff(kf_E)));
