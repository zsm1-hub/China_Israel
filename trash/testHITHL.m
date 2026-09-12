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
num_to_select = 7000;
% nparticles=10000; % numbers of particles
% num_to_select = 10000;
days=89.5;  % days
seconds=20.0;  % days
dt=0.1; % s  Advection_RK4 delta_t drift时间间隔
timerange=50:150;
Ttot = 0.1*(timerange(end)-timerange(1));
ini='box'
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
%% clear
p1=pairs_time;
clear pairs_time;
pairs_time=clear_nan_in_pairs_time(p1);
clear p1
%% error
% Ttot = days*24*3600;
Ttot = 100.*0.1;

tpts = length(pairs_time);
%
npairs = zeros(tpts,1);
for i = 1:tpts
    npairs(i) = length(find(~isnan(pairs_time(i).dul)));
end

%% Make into a single vector
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
    
end

%%
clear pairs_time

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



% Generate vel axis
vel_bins = linspace(-2, 2, 50);
vel_axis = 0.5*(vel_bins(1:end-1) + vel_bins(2:end));


%%
tic

pairs_sep = struct('dul', 'dut');

for i = 1:length(dist_axis)
    disp(i)
    id = find(dist>= dist_bin(i) & dist<dist_bin(i+1));
    
    pairs_sep(i).dul = dul(id);
    pairs_sep(i).dut = dut(id);
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

for i = 1:length(dist_axis)
    %pairs_per_bin(i) = length(id);
    %SF1l(i) = nanmean(pairs_sep(i).dul.^1);
    %SF1t(i) = nanmean(pairs_sep(i).dut.^1);
    
    SF2ll(i) = nanmean(pairs_sep(i).dul.^2);
    SF2tt(i) = nanmean(pairs_sep(i).dut.^2);
    %SF2lt(i) = nanmean(pairs_sep(i).dut.*pairs_sep(i).dul);
    
    %SF3lll(i) = nanmean(pairs_sep(i).dul.^3);
    %SF3ltt(i) = nanmean(pairs_sep(i).dul.*pairs_sep(i).dut.^2);
end
SF2=SF2ll+SF2tt;

%% test for setting up block bootstrap

test_flag =0 ;
if test_flag == 1
    ts = 1:12;
    blockSize = 2;
    numBlocks = length(ts) / blockSize;           % must be integer
    %blocks = reshape(ts, [numBlocks,blockSize])  % reshape into non-overlapping blocks
    blocks = reshape(ts, [blockSize, numBlocks])';
    nSamples = 10;
    samples = bootstrp(1, @(x)x', blocks);
    % the funny x' thing happens because the data is being converted to a row vector
    
end
%% Degree of freedom using time of process and total length of experiment
%
%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
Tscale_tot = 1./(((SF2ll +SF2tt).^0.5)./dist_axis);
Tscale_ll = 1./(((SF2ll).^0.5)./dist_axis);
Tscale_tt = 1./(((SF2tt).^0.5)./dist_axis);

dof = ceil(Ttot./Tscale_tot); % this is essentially T_tot/T_scale(r)

%%
for i = 1:length(pairs_sep)
    npairs_sep(i) = length(pairs_sep(i).dul);
    n_blocks_sep(i) =  dof(i); % number of blocks at that separation (basically the dof)
    nsamps_per_block_sep(i) = ceil(npairs_sep(i)/ n_blocks_sep(i));
end


clear SF3 SF3_mean SF3_stderr

num_boot = 1000;
SF3 = zeros(length(dist_axis), num_boot);
SF3_mean = zeros(length(dist_axis),1);
SF3_stderr = zeros(length(dist_axis),1);
SF1l = zeros(length(dist_axis), num_boot);
SF2l = zeros(length(dist_axis), num_boot);
SF3l = zeros(length(dist_axis), num_boot);
%%
tic
for i = 1:length(dist_axis)
    
    
    disp(i)
    blocksize = nsamps_per_block_sep(i);
    %blocksize = 1;
    numblocks = floor(npairs_sep(i)/ blocksize);
    
    if npairs_sep(i)>10
        n = numblocks*blocksize;
        
        blocks_dul = reshape(pairs_sep(i).dul(1:n), [blocksize, numblocks])';
        blocks_dut = reshape(pairs_sep(i).dut(1:n), [blocksize, numblocks])';
        
        SF3_samp = blocks_dul.^3 + blocks_dul.*blocks_dut.^2;
        SF1l_samp = blocks_dul;
        SF2l_samp = blocks_dul.^2;
        SF3l_samp = blocks_dul.^3;
        
        
        SF3(i,:) = bootstrp(num_boot, @(x)mean(mean(x,2),1), SF3_samp);
        SF1l(i,:) = bootstrp(num_boot, @(x)mean(mean(x,2),1), SF1l_samp);
        SF2l(i,:) = bootstrp(num_boot, @(x)mean(mean(x,2),1), SF2l_samp);
        SF3l(i,:) = bootstrp(num_boot, @(x)mean(mean(x,2),1), SF3l_samp);
        % the double mean above first takes mean over the blocks, then averages
        % the different blocks.
    else
        SF3(i,:) = NaN;
        SF1l(i,:) = NaN;
        SF2l(i,:) = NaN;
        SF3l(i,:) = NaN;
    end
    % Mean and standard error of the estimates
    SF3_mean(i) = mean(SF3(i,:));
    SF3_stderr(i) = std(SF3(i,:)); % boot strap std err is the std of bs estimates
    
end
toc
outputname=[input_dir,Case,'_pars_P',num2str(nparticles),'T',num2str(timerange(end)),...
    ini,'bootstrapHL2.mat']
% [input_dir,'wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
save(outputname,'SF3','SF3_mean', 'SF3_stderr', 'dof',...
     'dist_axis', 'dist_bin','SF2','SF2ll','SF2tt','SF1l','SF2l','SF3l','Th_all','nsample')
