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
Case='nowave'; % wave
nparticles=289; % numbers of particles
days=89.5;  % days
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
% ini='_roughbox200g_500m'
ini='_roughsmall'
timerange=1:1940;
% timerange=1:720;
% timerange=1220:1940;
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

save([fname(1:end-3),ini,'traj.mat'],'traj');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 3. Calc pairs of particles' variables 
%                    pair_time 储存不同粒子对的变量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% follow Balwada's routine : "trajectories2binnedpairs_vectorized.m"
% need code "dist_rx.m","dist_ry.m","dist_geo.m","dist_du.m"

Htraj=traj.H;
%% Break up by time and separation bins
% took 58.9s for GLAD
% took 146s for LASER 
tic 
for i=1:length(traj.T_axis)
    disp(i)
    id = find(Htraj(i,:)<-500);

    
    if length(id)>1
        X = traj.trajmat_X(i,id)';
        Y = traj.trajmat_Y(i,id)';
        U = traj.trajmat_U(i,id)';
        V = traj.trajmat_V(i,id)';
        
        Xvec = [X, Y];
        
        pairs_time(i).dist = pdist(Xvec, @dist_geo);
        
        rx = pdist(Xvec, @dist_rx);
        ry = pdist(Xvec, @dist_ry);
        
        magr = sqrt(rx.^2 + ry.^2);
        
        rx = rx./magr; ry = ry./magr;
        theta = atan2(ry, rx);
        
        dux = pdist(U, @dist_du);
        duy = pdist(V, @dist_du);
        
        
        pairs_time(i).dul = dux.*rx + duy.*ry;
        pairs_time(i).dut = duy.*rx - dux.*ry;
        pairs_time(i).theta=theta;
    else
        pairs_time(i).dul = NaN;
        pairs_time(i).dut = NaN;
        pairs_time(i).dist = NaN;
        pairs_time(i).theta = NaN;
    end
end
toc

p1=pairs_time;
clear pairs_time;
pairs_time=clear_nan_new(p1);
clear p1

Ttot = length(timerange)*3600;

tpts = length(pairs_time);

npairs = zeros(tpts,1);
for i = 1:tpts
    npairs(i) = length(find(~isnan(pairs_time(i).dul)));
end

%% Make into a single vector
%
dul = zeros(sum(npairs),1);
dut = zeros(sum(npairs),1);
dist = zeros(sum(npairs),1);
theta = zeros(sum(npairs),1);
%
% estimate num of pairs

empty1 = 1;
for i = 1:tpts % time loop
    if npairs(i) == 0
        continue
    end
    
    if npairs(i) == 1
        dist(empty1) = pairs_time(i).dist;
        dul(empty1) = pairs_time(i).dul;
        dut(empty1) = pairs_time(i).dut;
        theta(empty1) = pairs_time(i).theta;
        empty1 = empty1+1;
    end
    
    if npairs(i) >1
        dist(empty1: empty1+npairs(i)-1) = pairs_time(i).dist;
        dul(empty1: empty1+npairs(i)-1) = pairs_time(i).dul;
        dut(empty1: empty1+npairs(i)-1) = pairs_time(i).dut;
        theta(empty1: empty1+npairs(i)-1) = pairs_time(i).theta;
        empty1 = empty1+npairs(i);
    end
    
end

%%
clear pairs_time

if strcmpi(ini, '_roughsmall_500m') || strcmpi(ini, '_rough_500m') || strcmpi(ini, '_roughbox200g_500m') || strcmpi(ini, '_roughbox100g_500m')

    
    % gamma = 1.3;
    % dist_bin(1) = 1000; % in m
    gamma = 1.3;
    dist_bin(1) = 1000; % in m
    dist_bin = gamma.^[0:100]*dist_bin(1);
    id = find(dist_bin>600*10^3,1);
    dist_bin = dist_bin(1:id-1);
    dist_bin(2:end+1) = dist_bin(1:end);
    dist_bin(1) = 0;
    dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));
    
    
    % Generate vel axis
    vel_bins = linspace(-2, 2, 50);
    vel_axis = 0.5*(vel_bins(1:end-1) + vel_bins(2:end));
else
    gamma = 1.3;

    dist_bin(1) = 4000; % in m
    dist_bin = gamma.^[0:100]*dist_bin(1);
    id = find(dist_bin>600*10^3,1);
    dist_bin = dist_bin(1:id-1);
    dist_bin(2:end+1) = dist_bin(1:end);
    dist_bin(1) = 0;
    dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));
    
    
    % Generate vel axis
    vel_bins = linspace(-2, 2, 50);
    vel_axis = 0.5*(vel_bins(1:end-1) + vel_bins(2:end));
end

%%
tic

pairs_sep = struct('dul', 'dut','theta');

for i = 1:length(dist_axis)
    disp(i)
    id = find(dist>= dist_bin(i) & dist<dist_bin(i+1));
    
    pairs_sep(i).dul = dul(id);
    pairs_sep(i).dut = dut(id);
    pairs_sep(i).theta = theta(id);
end
toc
%% add by zsm
for ij=1:length(pairs_sep)
    nsample(ij)=length(pairs_sep(ij).dul);
end
save([Case,'P',num2str(nparticles),ini,'nsamp.mat'], ...
    'nsample','dist_axis','dist_bin');
%%
% Compute mean SF2 to use for estimating DOF
nbins = 16;
theta_edges = linspace(-pi, pi, nbins+1);
theta_mid   = theta_edges(1:end-1) + diff(theta_edges)/2;
dul_theta = NaN(length(dist_axis), nbins);
dut_theta = NaN(length(dist_axis), nbins);
for i = 1:length(dist_axis)
    
    dul_i   = pairs_sep(i).dul;
    dut_i   = pairs_sep(i).dut;
    theta_i = pairs_sep(i).theta;
    
    for j = 1:nbins
        idx = theta_i >= theta_edges(j) & ...
              theta_i <  theta_edges(j+1);
        
        if sum(idx) > 1
            dul_theta(i,j) = nanmean(dul_i(idx));
            dut_theta(i,j) = nanmean(dut_i(idx));
        else
            dul_theta(i,j) = NaN;
            dut_theta(i,j) = NaN;
        end
    end
end
outputname=[input_dir,Case,'_pars_P',num2str(nparticles),'T',num2str(timerange(end)),...
    ini,'theta.mat']
% [input_dir,'wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
save(outputname,'dul_theta','dut_theta',...
     'dist_axis', 'dist_bin','theta_mid')

%----------------------------------------------- end of test, change by zsm ---------------------------------


