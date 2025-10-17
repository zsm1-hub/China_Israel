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
nparticles=65536; % numbers of particles
num_to_select = 5000;
days=89.5;  % days
seconds=0.05;  % days
dt=2.5e-4; % s  Advection_RK4 delta_t drift时间间隔
timerange=1:180;
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
input_dir='/meddy/simingzhang/Data/Parcels_data/HIT2d_rough/';
addpath(input_dir)
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
time_batch=2;
num_workers = 4; % 例如使用4个工作节点
% HIT2d_pars_P2500T0.2seconds.nc
fname=[input_dir,'wave_pars_P',num2str(nparticles),'T',num2str(days),'seconds.nc'];
if strcmpi(Case, 'wave')
    fname=[input_dir,'wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
end

if strcmpi(Case, 'nowave')
    fname=[input_dir,'nowave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
end


if strcmpi(Case, 'HIT2d')
    fname=[input_dir,Case,'_pars_P',num2str(nparticles),'T',num2str(seconds),'seconds.nc'];
    oname=[input_dir,Case,'_pars_P',num2str(num_to_select),'T',num2str(seconds),'seconds.nc'];
end

lon=ncread(fname,'lon');
lat=ncread(fname,'lat');

ue=ncread(fname,'ue');
ve=ncread(fname,'ve');

num_columns = size(lon, 2); % 65536
random_indices = randperm(num_columns, num_to_select); % 随机选择不重复索引

% 抽取这些列
lons = lon(:, random_indices);
lats = lat(:, random_indices);
ues = ue(:, random_indices);
ves = ve(:, random_indices);

% read coarse-graining
xscale=[2:18,21:3:48,54:6:120,132:12:240,264:24:504];
PI=zeros(1,length(xscale));
for iii=1:length(xscale)
    pistr=['th',num2str(xscale(iii))];
    eval(['th',num2str(xscale(iii)),'=ncread(fname,','''',pistr,'''',');'])
    % eval(['Th(',num2str(iii),')=nanmean(','th',num2str(xscale(iii)),'(:));'])
    eval(['Th_all(',num2str(iii),',:)=nanmean(','th', ...
        num2str(xscale(iii)),'(1:183,','random_indices','),2);'])
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

lon=lon(timerange,:);
lat=lat(timerange,:);
u=ue(timerange,:);
v=ve(timerange,:);
%%%%%%%%%%%%%%%%%%%%%there is a 2D experiment, So I assuming H=-501
%%%%%%%%%%%%%%%%%%%%%H 没有意义, 只是在Balwada的code里面只采样了500米以上的粒子

traj=struct();
traj.trajmat_X=lon;traj.trajmat_Y=lat;
traj.trajmat_U=u;traj.trajmat_V=v;
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
for i=1:length(traj.T_axis)
    disp(i)
    %id = find(~isnan(traj.trajmat_X(i,:))); % find non-NaN
    
    %id = find(~isnan(traj.trajmat_X(i,:)) & Htraj(i,:)<-500); % find non-Nan and deep
    % id = find(~isnan(traj.trajmat_X(i,:)) & Htraj(i,:)<-500 & ...
    %     traj.trajmat_X(i,:)>=-91 & traj.trajmat_X(i,:)<=-84 & ...
    %     traj.trajmat_Y(i,:)>=24); % find non-Nan and deep and in similar region to GLAD
    id = find(Htraj(i,:)<-500);
    % if mod(i,300)==0
    %     disp(i)
    % end
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              4. Calc 2-order and 3-order structure function (time-mean)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% follow Balwada's routine : "pairtimes2SF.m"

tpts = length(pairs_time);
%
npairs = zeros(tpts,1);
for i = 1:tpts
    % modifide by zsm
    % npairs(i) = length(find(~isnan(pairs_time(i).dul)));
    npairs(i) = length(pairs_time(i).dist);
end

%% % Align pairs in a single vector
% 
dul = zeros(sum(npairs),1);
dut = zeros(sum(npairs),1);
dist = zeros(sum(npairs),1);

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
        empty1 = empty1+1;
    end
    
    if npairs(i) >1
        dist(empty1: empty1+npairs(i)-1) = pairs_time(i).dist;
        dul(empty1: empty1+npairs(i)-1) = pairs_time(i).dul;
        dut(empty1: empty1+npairs(i)-1) = pairs_time(i).dut;
        empty1 = empty1+npairs(i);
    end
    disp(i)
end

%%
clear pairs_time
% 
% gamma = 1.5;
% 
% dist_bin(1) = 10; % in m
% dist_bin = gamma.^[0:100]*dist_bin(1);
% 
% % dist_bin for cg
% % dist_bin=[1:18 21:3:48 54:6:114].*2e3;
% 
% id = find(dist_bin>1000*10^3,1);
% dist_bin = dist_bin(1:id);
% dist_bin(2:end+1) = dist_bin(1:end);
% dist_bin(1) = 0;
% dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));
% dist_axis=dist_axis(dist_axis>2e3);
% Generate vel axis
% vel_bins = linspace(-2, 2, 50);
% vel_axis = 0.5*(vel_bins(1:end-1) + vel_bins(2:end));
% min_resolution = 0.0123 * sqrt(2); % 最小分辨率 ≈ 0.0174
% max_scale = pi * sqrt(2);          % 最大尺度 π×√2 ≈ 4.44
% 
% % 计算分箱数量
% bin_width = min_resolution; % 固定分箱宽度
% num_bins = ceil(max_scale / bin_width); % 向上取整
% 
% % 创建分箱边界
% dist_bin = (0:num_bins) * bin_width;
% 
% % 确保不超过最大尺度
% if dist_bin(end) > max_scale
%     dist_bin(end) = max_scale;
% end
% 
% % 计算分箱中心
% dist_axis = 0.5 * (dist_bin(1:end-1) + dist_bin(2:end));
% num_bins = 20; % 总箱数
% log_min = log10(6.2832/512.*sqrt(2)); % 1mm
% log_max = log10(pi.*sqrt(2)); % 6.28m
% 
% dist_bin = logspace(log_min, log_max, num_bins);
% % dist_bin = [0, dist_bin]; % 添加0起点
% dist_axis = 0.5 * (dist_bin(1:end-1) + dist_bin(2:end));
% load HIT2d_Eul_r.mat
% dist_bin=dist_axis(1:1:end);
% clear dist_axis
% dist_axis = 0.5 * (dist_bin(1:end-1) + dist_bin(2:end));

% min_val = 1e-20 ;
% max_val = pi * sqrt(2);
% 
% % 设置点数（例如10个点）
% num_points = 100;
% 
% % 在对数空间生成等距点（使用自然对数）
% log_min = log(min_val);
% log_max = log(max_val);
% 
% % 生成对数空间等距点
% log_points = linspace(log_min, log_max, num_points);
% 
% % 通过指数还原为线性尺度
% arr = exp(log_points);
% dist_bin=arr;
% dist_axis = 0.5 * (dist_bin(1:end-1) + dist_bin(2:end));
% dr=dist_axis(2:end)-dist_axis(1:end-1);
% disp(dr)
% a=find(dr>0.0123);
% disp([num2str(a(1)),'    ',num2str(dist_axis(a(1)))])
% dist_bin=dist_bin(a(1):end);
% dist_axis=dist_axis(a(1):end);

% test linear 
arr = linspace(0.0123,pi.*sqrt(2),120)
dist_bin=arr;
dist_axis = 0.5 * (dist_bin(1:end-1) + dist_bin(2:end));


% min_val = 0.0123 * sqrt(2); % ≈0.0174
% max_val = pi * sqrt(2);     % ≈4.4429
% 
% % 设置总点数
% total_points = 30;
% 
% % 计算对数空间范围
% log_min = log(min_val);
% log_max = log(max_val);

% % 在小尺度区域增加点数密度
% % 使用分段对数分布：小尺度区域点数更多
% transition_point = 0.1; % 设置过渡点（可调整）
% log_transition = log(transition_point);
% 
% % 计算小尺度区域点数比例（占总点数的60%）
% small_scale_ratio = 0.2;
% num_small_points = round(total_points * small_scale_ratio);
% num_large_points = total_points - num_small_points;
% 
% % 生成小尺度区域的对数点（更密集）
% small_log_points = linspace(log_min, log_transition, num_small_points);
% 
% % 生成大尺度区域的对数点（较稀疏）
% large_log_points = linspace(log_transition, log_max, num_large_points);
% 
% % 合并点并确保唯一性
% all_log_points = unique([small_log_points, large_log_points]);
% 
% % 转换回线性空间
% dist_bin = exp(all_log_points);
% 
% % 计算轴点和间距
% dist_axis = 0.5 * (dist_bin(1:end-1) + dist_bin(2:end));
% dr = dist_axis(2:end) - dist_axis(1:end-1);

tic
for i = 1:length(dist_axis)
    disp(i)
    id = find(dist>= dist_bin(i) & dist<dist_bin(i+1));
    
    pairs_sep(i).dul = dul(id);
    pairs_sep(i).dut = dut(id);
    
    pairs_per_bin(i) = length(id);
    SF1l(i) = nanmean(pairs_sep(i).dul.^1);
    SF1t(i) = nanmean(pairs_sep(i).dut.^1);
        
    SF2ll(i) = nanmean(pairs_sep(i).dul.^2);
    SF2tt(i) = nanmean(pairs_sep(i).dut.^2);
    SF2lt(i) = nanmean(pairs_sep(i).dut.*pairs_sep(i).dul);
    
    SF3lll(i) = nanmean(pairs_sep(i).dul.^3);
    SF3ltt(i) = nanmean(pairs_sep(i).dul.*pairs_sep(i).dut.^2);
end
toc
SF3_mean=(SF3lll+SF3ltt)'

save([input_dir,'test.mat'],'SF3_mean','dist_axis')