

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
num_to_select = 10000;
% nparticles=10000; % numbers of particles
% num_to_select = 10000;
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

% read coarse-graining
xscale=[2:18,21:3:48,54:6:120,132:12:240,264:24:504];
PI=zeros(1,length(xscale));
for iii=1:length(xscale)
    pistr=['th',num2str(xscale(iii))];
    eval(['th',num2str(xscale(iii)),'=ncread(fname,','''',pistr,'''',');'])
    % eval(['Th(',num2str(iii),')=nanmean(','th',num2str(xscale(iii)),'(:));'])
    eval(['Th_all(',num2str(iii),',:)=nanmean(','th', ...
        num2str(xscale(iii)),'(timerange,','random_indices','),2);'])
end
%%%%%%%%%%%%%%%%check right?%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([oname(1:end-3),'traj.mat'],'lons','lats','ues','ves','Th_all')
% load HIT2d_pars_P289T0.05secondstraj.mat
% fname='s2sflux_spec_hit_tukey.0002.nc';
% ncdisp(fname)
% % Thm_Eulerian=ncread(fname,'Thm');
% filtscale=ncread(fname,'filtscale');
% filtscale=filtscale(1:end-1);
% figure(1)
% semilogx(filtscale,nanmean(Th_all,2))
% figure(2)
% dt=2.5e-4;
% for t=1:size(lons,1)-1
%     U(t,:)=(lons(t+1,:)-lons(t,:))./dt;
%     V(t,:)=(lats(t+1,:)-lats(t,:))./dt;
% end
% traj=150;
% plot(U(:,traj));hold on
% plot(0.5.*(ues(2:end,traj)+ues(1:end-1,traj)))
% plot(ues(traj,:))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear lon;clear lat;clear ue;clear ve;
lon=lons;lat=lats;ue=ues;ve=ves;
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

% lon=lon;
% lat=lat;
u=ue;
v=ve;
%%%%%%%%%%%%%%%%%%%%%there is a 2D experiment, So I assuming H=-501
%%%%%%%%%%%%%%%%%%%%%H 没有意义, 只是在Balwada的code里面只采样了500米以上的粒子

traj=struct();
traj.trajmat_X=lon;traj.trajmat_Y=lat;
traj.trajmat_U=ue;traj.trajmat_V=ve;
traj.H=-520.*ones(size(v,1),size(v,2));
traj.T_axis=linspace(dt, (size(v,1))*dt, size(v,1))./86400;

save([fname(1:end-3),'traj.mat'],'traj');

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

% gamma = 1.3;
% 
% dist_bin(1) = 0.0123*2; % in m
% 
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

gamma = 1.5;
dist_bin(1) = 0.0123; % in m

dist_bin = gamma.^[0:100]*dist_bin(1);
%%% modifed by zsm
% id = find(dist_bin>3.14*sqrt(2),1);
id = find(dist_bin>6.28,1);

dist_bin = dist_bin(1:id-1);
dist_bin(2:end+1) = dist_bin(1:end);
% modified by zsm
dist_bin(1) = dist_bin(3)-dist_bin(2); 
dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));

dist_axis(2:end)-dist_axis(1:end-1)



% Generate vel axis
vel_bins = linspace(-2, 2, 50);
vel_axis = 0.5*(vel_bins(1:end-1) + vel_bins(2:end));

% load HIT2d_Eul_r.mat
% dist_bin=dist_axis(1:1:end);
% clear dist_axis
% dist_axis = 0.5 * (dist_bin(1:end-1) + dist_bin(2:end));

% 将外层循环改为 parfor
if isempty(gcp('nocreate'))
    parpool('local', num_workers);
end

parfor ii = 1:num_chunks
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
            pairs_time(i).dist = pdist(Xvec, @dist_nogeo);
            
            rx = pdist(Xvec, @dist_norx);
            ry = pdist(Xvec, @dist_nory);
            
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
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %              4. Calc 2-order and 3-order structure function (time-mean)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tpts = length(pairs_time);
    
    % 预分配 npairs 数组
    npairs = zeros(tpts, 1);
    for i = 1:tpts
        if ~isempty(pairs_time(i).dist)
            npairs(i) = length(pairs_time(i).dist);
        end
    end
    
    %% % Align pairs in a single vector
    total_pairs = sum(npairs);
    dul = zeros(total_pairs, 1);
    dut = zeros(total_pairs, 1);
    dist = zeros(total_pairs, 1);
    
    idx = 1;
    for i = 1:tpts
        if npairs(i) > 0
            end_idx = idx + npairs(i) - 1;
            dist(idx:end_idx) = pairs_time(i).dist;
            dul(idx:end_idx) = pairs_time(i).dul;
            dut(idx:end_idx) = pairs_time(i).dut;
            idx = end_idx + 1;
        end
    end
    

    
    %% 计算结构函数
    % 预分配结果数组
    pairs_per_bin = zeros(length(dist_axis), 1);
    SF1l = zeros(length(dist_axis), 1);
    SF1t = zeros(length(dist_axis), 1);
    SF2ll = zeros(length(dist_axis), 1);
    SF2tt = zeros(length(dist_axis), 1);
    SF2lt = zeros(length(dist_axis), 1);
    SF3lll = zeros(length(dist_axis), 1);
    SF3ltt = zeros(length(dist_axis), 1);
    nvaild = zeros(length(dist_axis), 1);
    
    for i = 1:length(dist_axis)
        % 找出距离在[dist_bin(i), dist_bin(i+1)]之间的点对
        idx_bin = dist >= dist_bin(i) & dist < dist_bin(i+1);
        nvaild(i) = sum(idx_bin);
        
        if nvaild(i) > 0
            dul_bin = dul(idx_bin);
            dut_bin = dut(idx_bin);
            
            SF1l(i) = nansum(dul_bin);
            SF1t(i) = nansum(dut_bin);
            SF2ll(i) = nansum(dul_bin.^2);
            SF2tt(i) = nansum(dut_bin.^2);
            SF2lt(i) = nansum(dut_bin .* dul_bin);
            SF3lll(i) = nansum(dul_bin.^3);
            SF3ltt(i) = nansum(dul_bin .* dut_bin.^2);
        end
    end
    toc
    
    % 使用辅助函数保存结果
    saveResults(fname, ii, dist_axis, nvaild, ...
                SF1l, SF1t, SF2ll, SF2tt, ...
                SF3ltt, SF3lll, time_batch, total_steps);
end
delete(gcp('nocreate'));

%% join
np=nparticles;
for jj=1:length(np)
    npart=np(jj);
    % SF1l_time=0;
    % SF1t_time=0;
    % SF2ll_time=0;
    % SF2tt_time=0;
    % SF3ltt_time=0;
    % SF3lll_time=0;
    % nvaild_time=0;

    % 集合之前分开计算的chunk
    ii=1;
    eval(['load ',input_dir,Case,'_pars_P',num2str(npart),'T20.0secondschunk',num2str(ii,'%02d'),'SF123.mat']);
    num_chunk=floor(total_steps./time_batch);
    remainder=mod(total_steps,time_batch);

    for ii = 1:num_chunk
        eval(['load ',input_dir,Case,'_pars_P',num2str(npart),'T20.0secondschunk',num2str(ii,'%02d'),'SF123.mat'])

        SF1t_time(:,ii)=SF1t./nvaild;
        SF1l_time(:,ii)=SF1l./nvaild;
        SF2tt_time(:,ii)=SF2tt./nvaild;
        SF2ll_time(:,ii)=SF2ll./nvaild;
        SF3ltt_time(:,ii)=SF3ltt./nvaild;
        SF3lll_time(:,ii)=SF3lll./nvaild;
        nvaild_time(:,ii)=nvaild;
        disp(ii);
    end

  

    % 计算平均值
    SF1t_time(SF1t_time==0)=nan;
    SF1l_time(SF1l_time==0)=nan;


    SF2tt_time(SF2tt_time==0)=nan;
    SF2ll_time(SF2ll_time==0)=nan;

    SF3ltt_time(SF3ltt_time==0)=nan;
    SF3lll_time(SF3lll_time==0)=nan;

end

SF3_mean=nanmean((SF3lll_time+SF3ltt_time),2);
SF3_mean./dist_axis'
dist_axis
% outputname=[input_dir,Case,'_pars_P',num2str(nparticles),'T',num2str(timerange(end)),...
%     'timeavg2.mat']
% outputname=[input_dir,Case,'_pars_P',num2str(nparticles),'T',num2str(timerange(end)),...
%     'timeavg3.mat']
outputname=[input_dir,Case,'_pars_P',num2str(nparticles),'T',num2str(timerange(end)),...
    'timeavg7.mat']
save(outputname,'SF3_mean','dist_axis','Th_all','nvaild_time')

dist_axis(2:end)-dist_axis(1:end-1)

clear;
% load test2500_1.mat
% load test.mat
% load HIT2d_pars_P10000T150timeavg2.mat
% load HIT2d_pars_P22500T150timeavg2.mat
% load HIT2d_pars_P62500T150timeavg2.mat
% load HIT2d_pars_P10000T150timeavg3.mat
% load HIT2d_pars_P22500T150timeavg3.mat
% load HIT2d_pars_P10000T150timeavg4.mat
% load HIT2d_pars_P22500T150timeavg4.mat
% load HIT2d_pars_P10000T150timeavg5.mat
load HIT2d_pars_P62500T150timeavg8.mat


kf1=1./dist_axis'.*2.*pi

% version 3
% str1=0.05;
% en1=5;
% en1=2.7;
% lambda_opt_idx=[100,10,1,1e-1];

% str1=0.02;
% en1=2.7;

% for version 4
% str1=0.05;
% en1=5;
% lambda_opt_idx=[100,10,1,1e-1];


% for version 4 for 62500 (30000)
str1=0.012;
en1=5;
% lambda_opt_idx=[1000,100,10,1,1e-1,1e-2];
% lambda_opt_idx=[100,10,1,1e-1];
% lambda_opt_idx=[1000,100,10,1,1e-1,1e-2];
lambda_opt_idx=[1e-5];



% str1=0.01;
% en1=5;

figure(3)
semilogx(dist_axis,SF3_mean./dist_axis')
ylim([0,20])

point=1
SF3_coarse = movmean(SF3_mean, point);
dist_axis_coarse=movmean(dist_axis, point);

figure(3)
semilogx(dist_axis_coarse,SF3_coarse./dist_axis_coarse')
ylim([0,20])
SF3_mean=SF3_coarse;
dist_axis=dist_axis_coarse;
% [SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(SF3_mean,dist_axis,0.009,6.32, ...
%     'log','RLS',1e-10);
% [residual_norms, solution_norms,...
% lambda_opt_idx]=Fk_fitting_SF3_Lcurve(SF3_mean(1:end),dist_axis(1:end),str1,en1, ...
%     'fuc','RLS',[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
%     1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14,1e-15],kf1);

% lambda_opt_idx=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
%     1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14,1e-15];
% clf
kf1= logspace(log10(1/5), log10(1/0.04), 9).*2.*pi;

[SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3_Lcurve(SF3_mean(1:end),dist_axis(1:end),str1,en1, ...
    'fuc','RLS',lambda_opt_idx,kf1,0);

figure(2)
semilogx(dist_axis,SF3_mean./dist_axis');hold on
semilogx(lf,Vt./lf');
% ylim([0,20])

figure(1)
semilogx(1./kf.*2.*pi,SpecFlux)
% ylim([-12,5])

% 

% load HIT2d_pars_P10000T89.5bootstrap.mat
% figure(2)
% semilogx(dist_axis,SF3_mean./dist_axis');hold on

% figure(4)
% loglog(dist_axis,SF3_mean);hold on
