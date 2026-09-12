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
Case='nowave'; % wave
nparticles=625; % numbers of particles
days=89.5;  % days
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
ini='_rough'
if strcmpi(ini, '_grid')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
end
if strcmpi(ini, '_rough')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_roughdistr_tukey/';
end
% input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
timerange=1:2140;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              4. Calc 2-order and 3-order structure function (time-mean)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% follow Balwada's routine : "pairtimes2SF.m"

tpts = length(pairs_time);
%
npairs = zeros(tpts,1);
for i = 1:tpts
    % npairs(i) = length(find(~isnan(pairs_time(i).dul)));
    npairs(i) = length(pairs_time(i).dist);;
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

gamma = 1.5;

dist_bin(1) = 10; % in m
dist_bin = gamma.^[0:100]*dist_bin(1);

% dist_bin for cg
% dist_bin=[1:18 21:3:48 54:6:114].*2e3;

id = find(dist_bin>1000*10^3,1);
dist_bin = dist_bin(1:id);
dist_bin(2:end+1) = dist_bin(1:end);
dist_bin(1) = 0;
dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));
% dist_axis=dist_axis(dist_axis>2e3);
% Generate vel axis
% vel_bins = linspace(-2, 2, 50);
% vel_axis = 0.5*(vel_bins(1:end-1) + vel_bins(2:end));


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              5. plot 2-order and 3-order structure function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Some exploratory plots below
% SF2
% 
% figure
% loglog(dist_axis, SF2ll, 'linewidth',2)
% hold all
% loglog(dist_axis, SF2tt, 'linewidth',2)
% 
% loglog(dist_axis, SF2ll+SF2tt, 'linewidth',2)
% 
% 
% loglog(dist_axis, 1e-4*dist_axis.^(2/3), '--', 'color','k')
% loglog(dist_axis, 1e-7*dist_axis.^(2), '--', 'color','k')

%axis([10 1000e3 5e-5 5e-1])

%% SF3 
% figure
% loglog(dist_axis, abs(SF3lll+SF3ltt), 'linewidth',2)
% hold all
% loglog(dist_axis, SF3lll+SF3ltt, '+', 'linewidth',2)
% loglog(dist_axis, -SF3lll-SF3ltt, 'o', 'linewidth',2)
% axis([10 1000e3 1e-8 1])

% save([fname(1:end-3),'SF.mat'],'dist_axis','SF2ll','SF2tt','SF3ltt','SF3lll');
% save([fname(1:end-3),'wholetime','SF.mat'],'dist_axis','SF2ll','SF2tt','SF3ltt','SF3lll');

% save([fname(1:end-3),'SForigin.mat'],'dist_axis','pairs_sep', ...
%     '-v7.3');
save([fname(1:end-3),'SF123',ini,'.mat'],'dist_axis', ...
    'SF1l','SF1t','SF2ll','SF2tt','SF3ltt','SF3lll','-v7.3');
% save([fname(1:end-3),'SF123_c.mat'],'dist_axis', ...
%     'SF1l','SF1t','SF2ll','SF2tt','SF3ltt','SF3lll','-v7.3');
%% 
load wave_pars_P289T89.5daysSF123.mat
loglog(dist_axis, abs(SF3lll+SF3ltt), 'linewidth',2)
hold all
loglog(dist_axis, SF3lll+SF3ltt, '+', 'linewidth',2)
loglog(dist_axis, -SF3lll-SF3ltt, 'o', 'linewidth',2)
load wave_pars_P625T89.5daysSF123.mat
loglog(dist_axis, abs(SF3lll+SF3ltt), 'linewidth',2)
hold all
loglog(dist_axis, SF3lll+SF3ltt, '+', 'linewidth',2)
loglog(dist_axis, -SF3lll-SF3ltt, 'o', 'linewidth',2)