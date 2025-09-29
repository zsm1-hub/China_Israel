%   data source: Iceland (wave and no wave case)
%   utility: import parcels data to calc SF2,SF3
%   doesn't use bootstrap to resample,insteadly, calc time-mean SF2 and SF3 directly
%   code writer: zsm, modified from Balwada 2022 sciadv supplyment
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
nparticles=15376; % numbers of particles
days=89.5;  % days
ini='_grid'; % grid rough repeat
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
% input_dir='/meddy/simingzhang/Data/Parcels_data/onetime_spectukey/';
if strcmpi(ini, '_grid')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
end
if strcmpi(ini, '_rough')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_roughdistr_tukey/';
end

addpath(input_dir)
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
time_batch=1;

if strcmpi(Case, 'wave')
    fname=[input_dir,'wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
end

if strcmpi(Case, 'nowave')
    fname=[input_dir,'nowave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
end

lon=ncread(fname,'lon');
lat=ncread(fname,'lat');

ue=ncread(fname,'ue');
ve=ncread(fname,'ve');


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

% save([fname(1:end-3),'traj.mat'],'traj');

%%%%%%%%%%%%%%%%%%% 时间维度按time_batchsize分块，防止内存溢出
TRAJ=traj;
total_steps = size(traj.trajmat_X, 1);  % 511
num_chunks = floor(total_steps / time_batch);  % 10
remainder = mod(total_steps, time_batch);      % 11
clear traj;


for i = 1:num_chunks
    % 当前分块的索引范围
    start_idx = (i-1)*time_batch + 1;
    end_idx = i*time_batch;
    
    % 提取当前分块的数据
    traj = struct();
    traj.trajmat_X = TRAJ.trajmat_X(start_idx:end_idx, :);
    traj.trajmat_Y = TRAJ.trajmat_Y(start_idx:end_idx, :);
    traj.trajmat_U = TRAJ.trajmat_U(start_idx:end_idx, :);
    traj.trajmat_V = TRAJ.trajmat_V(start_idx:end_idx, :);
    traj.H = TRAJ.H(start_idx:end_idx, :);
    traj.T_axis = TRAJ.T_axis(start_idx:end_idx);
    
    % 保存为MAT文件
    % save([fname(1:end-3),'traj.mat'],'traj');
    save([fname(1:end-3),'chunk',num2str(i,'%02d'),'traj.mat'], 'traj');
    fprintf('Saved chunk %d (steps %d-%d)\n', i, start_idx, end_idx);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 3. Calc pairs of particles' variables 
%                    pair_time 储存不同粒子对的变量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% follow Balwada's routine : "trajectories2binnedpairs_vectorized.m"
% need code "dist_rx.m","dist_ry.m","dist_geo.m","dist_du.m"



gamma = 1.5;
dist_bin = gamma.^(0:100) * 10; % 初始距离箱
id_stop = find(dist_bin > 1000e3, 1); % 找到第一个超过1000km的箱子
if ~isempty(id_stop)
    dist_bin = dist_bin(1:id_stop);
end
dist_bin = [0, dist_bin]; % 在开头插入0
dist_axis = 0.5 * (dist_bin(1:end-1) + dist_bin(2:end));

num_workers = 4; % 例如使用4个工作节点

if isempty(gcp('nocreate'))
    parpool('local', num_workers);
end

% 将外层循环改为 parfor
parfor ii = 1:num_chunks+1
    % 使用函数加载数据以避免 eval 在并行环境中的问题
    traj = loadTrajData(fname, ii);
    
    Htraj = traj.H;
    tic 
    
    % 预分配 pairs_time 数组
    pairs_time = struct('dist', cell(1, length(traj.T_axis)), ...
                       'dul', cell(1, length(traj.T_axis)), ...
                       'dut', cell(1, length(traj.T_axis)));
    
    for i = 1:length(traj.T_axis)
        % 使用逻辑索引替代 find 提高效率
        id = Htraj(i,:) < -500;
        
        if sum(id) > 1
            X = traj.trajmat_X(i,id)';
            Y = traj.trajmat_Y(i,id)';
            U = traj.trajmat_U(i,id)';
            V = traj.trajmat_V(i,id)';
            
            Xvec = [X, Y];
            
            % 计算距离
            pairs_time(i).dist = pdist(Xvec, @dist_geo);
            
            rx = pdist(Xvec, @dist_rx);
            ry = pdist(Xvec, @dist_ry);
            
            magr = sqrt(rx.^2 + ry.^2);
            
            rx = rx./magr; 
            ry = ry./magr;
            
            dux = pdist(U, @dist_du);
            duy = pdist(V, @dist_du);
            
            pairs_time(i).dul = dux.*rx + duy.*ry;
            pairs_time(i).dut = duy.*rx - dux.*ry;
        else
            pairs_time(i).dul = NaN;
            pairs_time(i).dut = NaN;
            pairs_time(i).dist = NaN;
        end
    end
    toc
    saveResults_pairs(fname, ii, ...
                pairs_time);
end
