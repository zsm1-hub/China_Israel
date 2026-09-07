%   data source: Iceland (wave and no wave case)
%   utility: import parcels data to calc SF2,SF3
%   doesn't use bootstrap to resample,insteadly, calc time-mean SF2 and SF3 directly
%   code writer: zsm, modified from Balwada 2022 sciadv supplyment

%%%%%%%%%%%% test dist_bin code ########################
clear all;close all;clc
% addpath('D:\LIN2023\model\RoyBarkan\LLC4320/')
% addpath('D:\LIN2023\crocotools\Preprocessingtools') % add function "spheric_dist.m"
% 
addpath('/meddy/simingzhang/Analysis/matlab/Parcels_SF/')
addpath('/meddy/simingzhang/Data/Parcels_data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          1. Basic setup and read data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Case='wave'; % wave
nparticles=289; % numbers of particles
days=89.5;  % days
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
ini='_roughsmall'
timerange=1:1940;

nblock_I = 4;
nblock_J = 4;
min_pairs = 500;


%%% useless
if strcmpi(ini, '_grid')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
    xscale=[2:18,21:3:48,54:6:114];
end
%%% 2km whole grid
if strcmpi(ini, '_rough') || strcmpi(ini, '_rough_1mon') || strcmpi(ini, '_rough_2mon')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_roughdistr_tukey/';
    xscale=[2:18,21:3:48,54:6:114];
end
%%% 2km 70*70box~140km
if strcmpi(ini, '_roughsmall') || strcmpi(ini, '_roughsmall_1mon') || strcmpi(ini, '_roughsmall_2mon') || strcmpi(ini, '_roughsmall_3mon')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_roughsmallregion/';
    xscale=[2:18,21:3:48,54:6:114];
end
%%% 2km LASER but in smallregion
if strcmpi(ini, '_roughLASER')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_roughLASER/';
    xscale=[2:18,21:3:48,54:6:114];
end
%%% 2month 140km box but in smallregion
% if strcmpi(ini, '_2month_roughsmall')
%     input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_2month_roughsmall/';
%     xscale=[2:18,21:3:48,54:6:114];
% end

%%% 2km 70*70box~140km rot
if strcmpi(ini, '_roughsmall_rot')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_roughsmallregion_rot/';
    xscale=[2:18,21:3:48,54:6:114];
end

if strcmpi(ini, '_roughsmall_div')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_roughsmallregion_div/';
    xscale=[2:18,21:3:48,54:6:114];
end
%%% 2km 70*70box~140km
if strcmpi(ini, '_cruise')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_cruise_roughsmallregion/';
    xscale=[2:18,21:3:48,54:6:114];
end
%%% 500 m 280*280box~140km
if strcmpi(ini, '_roughsmall_500m') 
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_roughsmallregion_500m/';
    xscale=[2:18,21:3:48,54:6:114,120:12:228,240:24:336];
end
%%% 500 m whole grid
if strcmpi(ini, '_rough_500m')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_rough_500m/';
    xscale=[2:18,21:3:48,54:6:114,120:12:228,240:24:336];
end
%%% 500 m 200*200box~100km
if strcmpi(ini, '_roughbox200g_500m')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_roughbox200g_500m/';
    xscale=[2:18,21:3:48,54:6:114,120:12:228,240:24:336];
end
%%% 500 m 100*100box~50km
if strcmpi(ini, '_roughbox100g_500m')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_roughbox100g_500m/';
    xscale=[2:18,21:3:48,54:6:114,120:12:228,240:24:336];
end
% input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
% timerange=1:2140;
% timerange=1:960;
% timerange=1:1200;
% timerange=1:720;
% timerange=1:1428;
% timerange=1:720;


if strcmpi(Case, 'wave')
    fname=[input_dir,'wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
end

if strcmpi(Case, 'nowave')
    fname=[input_dir,'nowave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
end

lon=ncread(fname,'lon');
lat=ncread(fname,'lat');

% ue=ncread(fname,'ue').*1852.*60.*cos(lat.*pi./180);
% ve=ncread(fname,'ve').*1852.*60;

ue=ncread(fname,'ue');
ve=ncread(fname,'ve');

lon=lon(timerange,:);
lat=lat(timerange,:);
ue=ue(timerange,:);
ve=ve(timerange,:);

% xscale=[2:18,21:3:48,54:6:114];
PI=zeros(1,length(xscale));
for iii=1:length(xscale)
    pistr=['th',num2str(xscale(iii))];
    eval(['th',num2str(xscale(iii)),'=ncread(fname,','''',pistr,'''',');'])
    % eval(['Th(',num2str(iii),')=nanmean(','th',num2str(xscale(iii)),'(:));'])
    eval(['Th_all(',num2str(iii),',:)=nanmean(','th', ...
        num2str(xscale(iii)),'(timerange,:),2);'])
end
%%%%%%%%%%%%%%%%check right?%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save([oname(1:end-3),'traj.mat'],'lons','lats','ues','ves','Th_all')

% read coarse-graining
% xscale=[2,4,6,8,10,12,16,20,30,50,60,100];
% PI=zeros(1,length(xscale));
% for iii=1:length(xscale)
%     pistr=['pi',num2str(xscale(iii))];
%     eval(['pi',num2str(xscale(iii)),'=ncread(fname,','''',pistr,'''',');'])
%     eval(['PI(',num2str(iii),')=nanmean(','pi',num2str(xscale(iii)),'(:));'])
% end
% semilogx(xscale,PI)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              2. Calc Lagrangian Velocity and save *traj.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for t=1:size(lon,1)-1
    U(t,:)=(spheric_dist(lat(t,:),lat(t,:),lon(t,:),lon(t+1,:)))./dt.*...
        sign(lon(t+1,:)-lon(t,:));
    V(t,:)=(spheric_dist(lat(t+1,:),lat(t,:),lon(t,:),lon(t,:)))./dt.*...
        sign(lat(t+1,:)-lat(t,:));
end
lon(end,:)=[];lat(end,:)=[];

lon=lon;
lat=lat;
u=U;
v=V;
%%%%%%%%%%%%%%%%%%%%%there is a 2D experiment, So I assuming H=-501
%%%%%%%%%%%%%%%%%%%%%H 没有意义, 只是在Balwada的code里面只采样了500米以上的粒子

traj=struct();
traj.trajmat_X=lon;traj.trajmat_Y=lat;
traj.trajmat_U=u;traj.trajmat_V=v;
traj.H=-520.*ones(size(v,1),size(v,2));
traj.T_axis=linspace(dt, (size(v,1))*dt, size(v,1))./86400;
%%%%%%%%%%%%%%%%%%%% read grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/meddy/simingzhang/Data/RB_iceland_data'))
gname='/meddy/simingzhang/Data/RB_iceland_data/niskin2km_500m_grd.nc'
lon_rho=ncread(gname,'lon_rho');
lat_rho=ncread(gname,'lat_rho');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nI,nJ] = size(lon_rho);

%% Build lon/lat -> continuous I/J mapping

[Igrid,Jgrid] = ndgrid(1:nI,1:nJ);

valid_grid = ...
    isfinite(lon_rho) & ...
    isfinite(lat_rho);

FI = scatteredInterpolant( ...
    lon_rho(valid_grid), ...
    lat_rho(valid_grid), ...
    Igrid(valid_grid), ...
    'linear', ...
    'nearest');

FJ = scatteredInterpolant( ...
    lon_rho(valid_grid), ...
    lat_rho(valid_grid), ...
    Jgrid(valid_grid), ...
    'linear', ...
    'nearest');

%% Check mapping using initial particle positions

I_initial = FI(lon(1,:),lat(1,:));
J_initial = FJ(lon(1,:),lat(1,:));


%% Define blocks in grid-index space

I_edges = round( ...
    linspace(1,nI+1,nblock_I+1));

J_edges = round( ...
    linspace(1,nJ+1,nblock_J+1));

nBlock = nblock_I*nblock_J;

fprintf('Grid: %d x %d\n',nI,nJ);
fprintf('Blocks: %d x %d\n',nblock_I,nblock_J);
fprintf('I block edges:\n');
disp(I_edges);
fprintf('J block edges:\n');
disp(J_edges);

%% Use your existing distance bins if available

if ~exist('dist_bin','var')

    gamma = 1.3;
    rmin = 4000;
    rmax = 600e3;

    dist_bin = rmin * gamma.^(0:100);
    dist_bin = dist_bin(dist_bin <= rmax);
    dist_bin = unique([0,dist_bin]);
end

dist_axis = 0.5 * ...
    (dist_bin(1:end-1) + dist_bin(2:end));

nScale = length(dist_axis);

%% Pair indices generated only once

nParticle = size(lon,2);

[pair_i,pair_j] = ...
    find(triu(true(nParticle),1));

%% Streaming accumulators

sum_SF1 = zeros(nScale,nBlock);
sum_SF2 = zeros(nScale,nBlock);
sum_SF3 = zeros(nScale,nBlock);

pair_count = zeros(nScale,nBlock);

%% Main loop

nTime = size(lon,1);

for it = 1:nTime

    if mod(it,100) == 0
        fprintf('time %d/%d\n',it,nTime);
    end

    valid_particle = ...
        isfinite(lon(it,:)) & ...
        isfinite(lat(it,:)) & ...
        isfinite(U(it,:)) & ...
        isfinite(V(it,:));

    valid_pair = ...
        valid_particle(pair_i) & ...
        valid_particle(pair_j);

    if ~any(valid_pair)
        continue;
    end

    pi = pair_i(valid_pair);
    pj = pair_j(valid_pair);

    lon1 = lon(it,pi)';
    lon2 = lon(it,pj)';

    lat1 = lat(it,pi)';
    lat2 = lat(it,pj)';

    u1 = U(it,pi)';
    u2 = U(it,pj)';

    v1 = V(it,pi)';
    v2 = V(it,pj)';

    %% Pair distance

    d = spheric_dist( ...
        lat1,lat2,lon1,lon2);

    %% Pair direction

    rx = spheric_dist( ...
        lat1,lat1,lon1,lon2);

    ry = spheric_dist( ...
        lat1,lat2,lon1,lon1);

    magr = sqrt(rx.^2 + ry.^2);

    good = ...
        isfinite(d) & ...
        isfinite(magr) & ...
        magr > 0;

    if ~any(good)
        continue;
    end

    d = d(good);

    rx = rx(good) ./ magr(good);
    ry = ry(good) ./ magr(good);

    %% Velocity increments

    dux = u2(good) - u1(good);
    duy = v2(good) - v1(good);

    dul = dux .* rx + duy .* ry;
    dut = duy .* rx - dux .* ry;

    %% Pair midpoint

    lon_mid = 0.5 * ...
        (lon1(good) + lon2(good));

    lat_mid = 0.5 * ...
        (lat1(good) + lat2(good));

    %% Convert midpoint to grid-index coordinates

    I_mid = FI(lon_mid,lat_mid);
    J_mid = FJ(lon_mid,lat_mid);

    %% Assign scale and spatial block

    scale_id = discretize(d,dist_bin);

    block_I = discretize(I_mid,I_edges);
    block_J = discretize(J_mid,J_edges);

    block_I(I_mid == I_edges(end)) = nblock_I;
    block_J(J_mid == J_edges(end)) = nblock_J;

    block_id = ...
        (block_I - 1) * nblock_J + block_J;

    good2 = ...
        isfinite(scale_id) & ...
        isfinite(block_id) & ...
        scale_id >= 1 & ...
        scale_id <= nScale & ...
        block_id >= 1 & ...
        block_id <= nBlock & ...
        isfinite(dul);

    if ~any(good2)
        continue;
    end

    scale_id = scale_id(good2);
    block_id = block_id(good2);

    dul_use = dul(good2);

    linear_id = sub2ind( ...
        [nScale,nBlock], ...
        scale_id,block_id);

    %% Accumulate SF1

    tmp = accumarray( ...
        linear_id, ...
        dul_use, ...
        [nScale*nBlock,1], ...
        @sum,0);

    sum_SF1 = sum_SF1 + ...
        reshape(tmp,nScale,nBlock);

    %% Accumulate SF2

    tmp = accumarray( ...
        linear_id, ...
        dul_use.^2, ...
        [nScale*nBlock,1], ...
        @sum,0);

    sum_SF2 = sum_SF2 + ...
        reshape(tmp,nScale,nBlock);

    %% Accumulate SF3

    % Longitudinal SF3:
    sf3_use = dul_use.^3;

    % If your paper uses the 2-D expression:
    %
    % dut_use = dut(good2);
    % sf3_use = dul_use.^3 + dul_use .* dut_use.^2;

    tmp = accumarray( ...
        linear_id, ...
        sf3_use, ...
        [nScale*nBlock,1], ...
        @sum,0);

    sum_SF3 = sum_SF3 + ...
        reshape(tmp,nScale,nBlock);

    %% Pair count

    tmp = accumarray( ...
        linear_id, ...
        1, ...
        [nScale*nBlock,1], ...
        @sum,0);

    pair_count = pair_count + ...
        reshape(tmp,nScale,nBlock);
end

%% Local structure functions

local_SF1 = nan(nScale,nBlock);
local_SF2 = nan(nScale,nBlock);
local_SF3 = nan(nScale,nBlock);

valid = pair_count >= min_pairs;

local_SF1(valid) = ...
    sum_SF1(valid) ./ pair_count(valid);

local_SF2(valid) = ...
    sum_SF2(valid) ./ pair_count(valid);

local_SF3(valid) = ...
    sum_SF3(valid) ./ pair_count(valid);

%% Calculate homogeneity

[H1,mean_SF1,std_SF1,rms_SF1,nvalid1] = ...
    calculate_H(local_SF1,pair_count,min_pairs);

[H2,mean_SF2,std_SF2,rms_SF2,nvalid2] = ...
    calculate_H(local_SF2,pair_count,min_pairs);

[H3,mean_SF3,std_SF3,rms_SF3,nvalid3] = ...
    calculate_H(local_SF3,pair_count,min_pairs);

outputname=[input_dir,Case,'_pars_P',num2str(nparticles),'T',num2str(timerange(end)),...
    ini,'gridBlockHL.mat']
% [input_dir,'wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
save(outputname,'H1','mean_SF1','std_SF1','rms_SF1','nvalid1',...
      'H2','mean_SF2','std_SF2','rms_SF2','nvalid2',...
      'H3','mean_SF3','std_SF3','rms_SF3','nvalid3',...
     'dist_axis')
