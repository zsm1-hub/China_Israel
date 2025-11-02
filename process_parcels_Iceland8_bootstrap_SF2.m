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
ini='_roughsmall_500m'
if strcmpi(ini, '_grid')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
    xscale=[2:18,21:3:48,54:6:114];
end
if strcmpi(ini, '_rough')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_roughdistr_tukey/';
    xscale=[2:18,21:3:48,54:6:114];
end
if strcmpi(ini, '_roughsmall')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_roughsmallregion/';
    xscale=[2:18,21:3:48,54:6:114];
end
if strcmpi(ini, '_roughsmall_rot')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_roughsmallregion_rot/';
    xscale=[2:18,21:3:48,54:6:114];
end

if strcmpi(ini, '_cruise')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_cruise_roughsmallregion/';
    xscale=[2:18,21:3:48,54:6:114];
end
if strcmpi(ini, '_roughsmall_500m')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_roughsmallregion_500m/';
    xscale=[1:18,21:3:48,54:6:114,120:12:228,240:24:336];
end

% input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
% timerange=1:2140;
% timerange=1:960;
% timerange=1:1200;
% timerange=1:720;
timerange=1:1940;

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
        
        pairs_time(i).dist = pdist(Xvec, @dist_geo);
        
        rx = pdist(Xvec, @dist_rx);
        ry = pdist(Xvec, @dist_ry);
        
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
Ttot = days*24*3600;
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
        
        SF2ll_samp = blocks_dul.^2; %+ blocks_dul.*blocks_dut.^2;
        SF2tt_samp = blocks_dut.^2;
        
        % create blocks of bootstrap samples
        %SF3_bs = bootstrp(num_boot, @(x)x', SF3_samp');
        % calculate means of each bootstrap sample
        %SF3 = mean(SF3_bs, 2);
        
        [SF2ll(i,:), bootsamp] = bootstrp(num_boot, @(x)mean(mean(x,2),1), SF2ll_samp);
        % the double mean above first takes mean over the blocks, then averages
        % the different blocks.
        % the below loop is needed because we want correspondence between
        % the SF2ll and SF2tt samples
        for j = 1:num_boot
            SF2tt(i,j) = mean(mean(SF2tt_samp(bootsamp(:,j),:),2),1);
        end
        
    else
        SF2ll(i,:) = NaN;
        SF2tt(i,:) = NaN; 
    end
    % Mean and standard error of the estimates
    SF2ll_mean(i)   = mean(SF2ll(i,:));
    SF2ll_stderr(i) = std(SF2ll(i,:)); % boot strap std err is the std of bs estimates
    SF2tt_mean(i)   = mean(SF2tt(i,:));
    SF2tt_stderr(i) = std(SF2tt(i,:));
end
toc
outputname=[input_dir,Case,'_pars_P',num2str(nparticles),'T',num2str(timerange(end)),...
    ini,'bootstrap_SF2.mat']
% [input_dir,'wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
save(outputname,'SF2ll','SF2tt', 'SF2ll_mean','SF2ll_stderr',...
    'dof','SF2tt_mean','SF2tt_stderr',...
     'dist_axis', 'dist_bin','nsample')

