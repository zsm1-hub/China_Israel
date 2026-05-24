clear all;close all;clc
% addpath('D:\LIN2023\model\RoyBarkan\LLC4320/')
% addpath('D:\LIN2023\crocotools\Preprocessingtools') % add function "spheric_dist.m"
% 
addpath('/meddy/simingzhang/Analysis/matlab/Parcels_SF/')
addpath('/meddy/simingzhang/Data/Parcels_data')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          1. Basic setup and read data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Case='HIT2d'; % wave
nparticles=62500; % numbers of particles
num_to_select = 1000;


days=89.5;  % days
seconds=20.0;  % days
dt=0.1; % s  Advection_RK4 delta_t drift时间间隔
timerange=50:150;
Ttot = 0.1*(timerange(end)-timerange(1));

% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
input_dir='/meddy/simingzhang/Data/Parcels_data/HIT2d_rough/';
addpath(input_dir)
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
time_batch=1;
num_workers = 2; % 例如使用4个工作节点
% HIT2d_pars_P2500T0.2seconds.nc
fname=[input_dir,'wave_pars_P',num2str(nparticles),'T',num2str(days),'seconds.nc'];
if strcmpi(Case, 'wave')
    fname=[input_dir,'wave_pars_P',num2str(nparticles),'T',num2str(days),'seconds.nc'];
end

if strcmpi(Case, 'nowave')
    fname=[input_dir,'nowave_pars_P',num2str(nparticles),'T',num2str(days),'seconds.nc'];
end


if strcmpi(Case, 'HIT2d')
    fname=[input_dir,Case,'_pars_P',num2str(nparticles),'T',num2str(seconds),'.0seconds.nc'];
    oname=[input_dir,Case,'_pars_P',num2str(num_to_select),'T',num2str(seconds),'seconds.nc'];
end

lon=ncread(fname,'lon');
lat=ncread(fname,'lat');

ue=ncread(fname,'ue');
ve=ncread(fname,'ve');

num_columns = size(lon, 2); % 65536
random_indices = randperm(num_columns, num_to_select); % 随机选择不重复索引

% 抽取这些列
lons = lon(timerange, random_indices);
lats = lat(timerange, random_indices);
ues = ue(timerange, random_indices);
ves = ve(timerange, random_indices);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              2. Calc Lagrangian Velocity and save *traj.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% for t=1:size(lon,1)-1
%     U(t,:)=(spheric_dist(lat(t,:),lat(t,:),lon(t,:),lon(t+1,:)))./dt.*...
%         sign(lon(t+1,:)-lon(t,:));
%     V(t,:)=(spheric_dist(lat(t+1,:),lat(t,:),lon(t,:),lon(t,:)))./dt.*...
%         sign(lat(t+1,:)-lat(t,:));
% end
% lon(end,:)=[];lat(end,:)=[];

lon=lons;
lat=lats;
ue=ues;
ve=ves;
%%%%%%%%%%%%%%%%%%%%%there is a 2D experiment, So I assuming H=-501
%%%%%%%%%%%%%%%%%%%%%H 没有意义, 只是在Balwada的code里面只采样了500米以上的粒子

traj=struct();
traj.trajmat_X=lon;traj.trajmat_Y=lat;
traj.trajmat_U=ue;traj.trajmat_V=ve;
traj.H=-520.*ones(size(ve,1),size(ve,2));
traj.T_axis=linspace(dt, (size(ve,1))*dt, size(ve,1))./86400;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 3. Calc pairs of particles' variables 
%                    pair_time 储存不同粒子对的变量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Htraj=traj.H;

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

        pairs_time(i).dist = pdist(Xvec, @dist_nogeo);
            
        rx = pdist(Xvec, @dist_norx);
        ry = pdist(Xvec, @dist_nory);
        
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 3. Calc pairs of particles' variables 
%                    pair_time 储存不同粒子对的变量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% follow Balwada's routine : "trajectories2binnedpairs_vectorized.m"
% need code "dist_rx.m","dist_ry.m","dist_geo.m","dist_du.m"

% gamma = 1.3;
% 
% dist_bin(1) = 0.0123*2; % in m

% gamma = 1.35;
% dist_bin(1) = 0.0123*sqrt(2)*2; % in m
% 
% dist_bin = gamma.^[0:100]*dist_bin(1);
% %%% modifed by zsm
% % id = find(dist_bin>3.14*sqrt(2),1);
% id = find(dist_bin>3.14,1);
% 
% dist_bin = dist_bin(1:id);
% dist_bin(2:end+1) = dist_bin(1:end);
% % modified by zsm
% dist_bin(1) = dist_bin(3)-dist_bin(2); 
% dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));

gamma = 1.2;
dist_bin(1) = 0.0123; % in m

dist_bin = gamma.^[0:100]*dist_bin(1);
%%% modifed by zsm
% id = find(dist_bin>3.14*sqrt(2),1);
id = find(dist_bin>3.14,1);

dist_bin = dist_bin(1:id-1);
dist_bin(2:end+1) = dist_bin(1:end);
% modified by zsm
dist_bin(1) = dist_bin(3)-dist_bin(2); 
dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));

nbins = 16;
theta_edges = linspace(-pi, pi, nbins+1);
theta_mid   = theta_edges(1:end-1) + diff(theta_edges)/2;


% Generate vel axis
vel_bins = linspace(-2, 2, 50);
vel_axis = 0.5*(vel_bins(1:end-1) + vel_bins(2:end));


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
            dull(i,j)=nanmean(dul_i(idx));
            dull2(i,j)=nanmean(dul_i(idx).^2);
            dull3(i,j)=nanmean(dul_i(idx).^3);
        else
            dull(i,j)=nan;
            dull2(i,j)=nan;
            dull3(i,j)=nan;
        end
    end
end

outputname=[input_dir,Case,'_pars_P',num2str(nparticles),'T',num2str(timerange(end)),...
    'theta.mat']
disp(outputname)
save(outputname,'dull','dull2','dull3','dist_axis','theta_mid')



