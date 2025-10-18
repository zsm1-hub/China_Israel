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
nparticles=289; % numbers of particles
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

lon=lon(timerange,:);
lat=lat(timerange,:);
ue=ue(timerange,:);
ve=ve(timerange,:);

xscale=[2:18,21:3:48,54:6:114];
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

gamma = 1.5;

dist_bin(1) = 10; % in m
dist_bin = gamma.^[0:100]*dist_bin(1);
id = find(dist_bin>1000*10^3,1);
dist_bin = dist_bin(1:id);
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
        
        SF3_samp = blocks_dul.^3 + blocks_dul.*blocks_dut.^2;
        
        % create blocks of bootstrap samples
        %SF3_bs = bootstrp(num_boot, @(x)x', SF3_samp');
        % calculate means of each bootstrap sample
        %SF3 = mean(SF3_bs, 2);
        
        SF3(i,:) = bootstrp(num_boot, @(x)mean(mean(x,2),1), SF3_samp);
        % the double mean above first takes mean over the blocks, then averages
        % the different blocks.
    else
        SF3(i,:) = NaN;
    end
    % Mean and standard error of the estimates
    SF3_mean(i) = mean(SF3(i,:));
    SF3_stderr(i) = std(SF3(i,:)); % boot strap std err is the std of bs estimates
    
end
toc
outputname=[input_dir,Case,'_pars_P',num2str(nparticles),'T',num2str(days),...
    'bootstrap.mat']
% [input_dir,'wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
save(outputname,'SF3','SF3_mean', 'SF3_stderr', 'dof',...
     'dist_axis', 'dist_bin','SF2','Th_all')

%% test plot
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
ini='_rough'
lambda=1e-10;

if strcmpi(Case, 'wave')
    cgname='s2sflux_spec_hf.0002.nc';
    cgstd='wavecase_modified_cg_tukey1_RodivofstrainofCG_Eul_std_.mat'
end

if strcmpi(Case, 'nowave')
    cgname='s2sflux_spec_smooth.0002.nc';
    cgstd='nowavecase_modified_cg_tukey1_RodivofstrainofCG_Eul_std_.mat'
end

% cgname='s2sflux_spec_smooth.0002.nc';
ncdisp(cgname)
Thm_Eulerian=ncread(cgname,'Thm');
filtscale=ncread(cgname,'filtscale');

outputname=[Case,'_pars_P',num2str(nparticles),'T',num2str(days),...
    'bootstrap.mat']

colors={'#0072BD','#D95319','#EDB120','#7E2F8E'};
colors_rgb = {...
    [0, 114/255, 189/255], ...   % #0072BD
    [217/255, 83/255, 25/255], ... % #D95319
    [237/255, 177/255, 32/255], ... % #EDB120
    [126/255, 47/255, 142/255] ... % #7E2F8E
};

load(outputname)

SF3=SF3;
% dstr=14;
dstr=1;
for tt=1:size(SF3,2)
     [SpecFlux(:,tt),Vt(:,tt),ebs(:,tt),kf,lf]=Fk_fitting_SF3(SF3(dstr:end,tt),dist_axis(dstr:end), ...
         2e3,200e3,'log','RLS',lambda);
end

clear CI_ebs CI_Vt CI_SpecFlux CI_SF3

% Convert parameters to variance preserving form for plotting
ebs_varp = ebs;
ebs_varp(1:end-1,:) = ebs(1:end-1,:).*kf';

for i = 1:size(ebs,1)
    CI_ebs(:,i) = prctile(ebs_varp(i,:), [95, 5]); 
end

for i = 1:size(Vt,1)
    CI_Vt(:,i) = prctile(Vt(i,:), [95,5]);
end

for i = 1:size(SF3,1)
    CI_SF3(:,i) = prctile(SF3(i,:), [95,5]);
end

for i = 1:size(SpecFlux,1)    
    CI_SpecFlux(:,i) = prctile(SpecFlux(i,:), [95,5]);
end

median_ebs = median(ebs_varp,2); 
mean_ebs = nanmean(ebs_varp,2); 
median_Vt = median(Vt,2);
mean_Vt = nanmean(Vt,2);

%median_S = median(S,2);
mean_SF3 = nanmean(SF3,2);

median_SpecFlux = median(SpecFlux,2);
mean_SpecFlux = nanmean(SpecFlux,2);


figure(1)
ii=2
nparticles=[289,625,2500,15376]; % numbers of particles

subplot(2,2,1)
semilogx(1./kf./1e3,nanmean(SpecFlux,2),'Marker','x','Color',colors{ii}, ...
    'LineWidth',1.5);
% semilogx(1./kf./1e3,(SpecFlux),'Marker','x','Color',colors{ii},'LineWidth',1.5);

hold on
% clear x_fill;clear y_fill
x_fill = [1./kf./1e3, fliplr(1./kf./1e3)];
y_fill = [CI_SpecFlux(1,:), fliplr(CI_SpecFlux(2,:))];
% y_fill = [(nanmean(SpecFlux,2)+std1)', fliplr((nanmean(SpecFlux,2)-std1)')];
fill(x_fill, y_fill, colors_rgb{ii}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

eval(['load ',Case,'_pars_P',num2str(nparticles(ii)),'T89.5daysCG_Lag_spectukey', ...
    ini,'.mat']);
std1=std(Th_all, 0, 2, 'omitnan');
x_fill = [filtscale'./1e3, fliplr(filtscale'./1e3)];
y_fill = [(Th'+std1)', fliplr((Th'-std1)')];
C{ii}=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{ii},'LineStyle','--');
fill(x_fill, y_fill, 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
errorbar(filtscale./1e3,Th', std1, 'Color', colors{ii})

grid on
ylim([-0.5e-7,0.5e-7]);
xlim([1e3,1e6]./1e3);
xlabel('km')
ylabel('m^{2}/s^{3}')
title([Case,': cg flux v.s. sf3-fitting fk'])
% legend([B{ii},C{ii}],{'SF3-fitting fk','Lag cg-spec fk'},'Location','southeast');
set(gca,'fontsize',16,'FontWeight','b')

subplot(2,2,2)
semilogx(dist_axis./1e3,nanmean(SF3,2)./dist_axis','Color',colors{ii}, ...
    'LineWidth',1.5);
hold on
semilogx(lf./1e3,nanmean(Vt,2)./lf','Marker','x','Color',colors{ii}, ...
    'LineWidth',1.5,'LineStyle','--');
grid on
% ylim([-4e-8,6e-8]);
xlim([1e3,1e6]./1e3);
xlabel('km')
ylabel('m^{2}/s^{3}')
title([Case,': sf3'])
% legend([B{ii},C{ii}],{'SF3-fitting fk','Lag cg-spec fk'},'Location','southeast');
set(gca,'fontsize',16,'FontWeight','b')


figure(2)
semilogx(1./kf./1e3,nanmean(ebs_varp(1:end-1,:),2),'Marker','x','Color',colors{ii}, ...
    'LineWidth',1.5);hold on
x_fill = [1./kf./1e3, fliplr(1./kf./1e3)];
y_fill = [CI_ebs(1,1:end-1), fliplr(CI_ebs(2,1:end-1))];
% y_fill = [(nanmean(SpecFlux,2)+std1)', fliplr((nanmean(SpecFlux,2)-std1)')];
fill(x_fill, y_fill, colors_rgb{ii}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
grid on
xlim([1e3,1e6]./1e3);
xlabel('km')
ylabel('m^{2}/s^{3}')
title([Case,': ebs*k'])
set(gca,'fontsize',16,'FontWeight','b')

% figure(1)
% subplot(2,2,3)
% fname{ii,:}=[Case,'_pars_P',num2str(nparticles(ii)),'T',num2str(days),'days.nc'];
% eval(['load ',fname{ii}(1:end-3),'SF123_alltime',ini,'.mat']);
% 
% SF3=(SF3lll_time+SF3ltt_time);
% dstr=14;
% dstr=1;
% for tt=1:size(SF3,2)
%      [SpecFlux(:,tt),Vt(:,tt),ebs(:,tt),kf,lf]=Fk_fitting_SF3(SF3(dstr:end,tt),dist_axis(dstr:end), ...
%          2e3,500e3,'log','RLS',lambda);
% end
% 
% std1=std(SpecFlux, 0, 2, 'omitnan');
% 
% semilogx(1./kf./1e3,nanmean(SpecFlux,2),'Marker','x','Color',colors{ii}, ...
%     'LineWidth',1.5);
% % semilogx(1./kf./1e3,(SpecFlux),'Marker','x','Color',colors{ii},'LineWidth',1.5);
% 
% hold on
% % clear x_fill;clear y_fill
% x_fill = [1./kf./1e3, fliplr(1./kf./1e3)];
% % y_fill = [CI_SpecFlux(1,:), fliplr(CI_SpecFlux(2,:))];
% y_fill = [(nanmean(SpecFlux,2)+std1)', fliplr((nanmean(SpecFlux,2)-std1)')];
% fill(x_fill, y_fill, colors_rgb{ii}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
% 
% grid on
% ylim([-0.5e-7,0.5e-7]);
% xlim([1e3,1e6]./1e3);
% xlabel('km')
% ylabel('m^{2}/s^{3}')
% title([Case,': cg flux v.s. sf3-fitting fk'])
% % legend([B{ii},C{ii}],{'SF3-fitting fk','Lag cg-spec fk'},'Location','southeast');
% set(gca,'fontsize',16,'FontWeight','b')
% 
% 
% subplot(2,2,4)
% semilogx(dist_axis./1e3,nanmean(SF3,2)./dist_axis','Color',colors{ii}, ...
%     'LineWidth',1.5);
% hold on
% semilogx(lf./1e3,nanmean(Vt,2)./lf','Marker','x','Color',colors{ii}, ...
%     'LineWidth',1.5,'LineStyle','--');
% grid on
% % ylim([-4e-8,6e-8]);
% xlim([1e3,1e6]./1e3);
% xlabel('km')
% ylabel('m^{2}/s^{3}')
% title([Case,': cg flux v.s. sf3-fitting fk'])
% % legend([B{ii},C{ii}],{'SF3-fitting fk','Lag cg-spec fk'},'Location','southeast');
% set(gca,'fontsize',16,'FontWeight','b')

figure(3)
semilogx(1./kf./1e3,nanmean(SpecFlux,2),'Marker','x','Color',colors{ii+1}, ...
    'LineWidth',1.5);hold on
load s625.mat
semilogx(1./kf./1e3,mean_SpecFlux,'Marker','x','Color',colors{ii}, ...
    'LineWidth',1.5);
load s289.mat
semilogx(1./kf./1e3,mean_SpecFlux,'Marker','x','Color',colors{ii-1}, ...
    'LineWidth',1.5);

grid on
xlim([1e3,1e6]./1e3);
xlabel('km')
ylabel('m^{2}/s^{3}')
title([Case,': fk'])
set(gca,'fontsize',16,'FontWeight','b')

%% plot 3
clear all;close all;clc
addpath('/meddy/simingzhang/Analysis/matlab/Parcels_SF/')
addpath('/meddy/simingzhang/Data/Parcels_data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          1. Basic setup and read data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Case='wave'; % wave
np=[289,625,2500]; % numbers of particles
days=89.5;  % days
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
ini='_rough'
lambda=1e-9;

if strcmpi(Case, 'wave')
    cgname='s2sflux_spec_hf_hann_corr.0002.nc';
    cgstd='wavecase_modified_cg_tukey1_RodivofstrainofCG_Eul_std_.mat'
end

if strcmpi(Case, 'nowave')
    cgname='s2sflux_spec_smooth.0002.nc';
    cgstd='nowavecase_modified_cg_tukey1_RodivofstrainofCG_Eul_std_.mat'
end

% cgname='s2sflux_spec_smooth.0002.nc';
ncdisp(cgname)
Thm_Eulerian=ncread(cgname,'Thm');
filtscale=ncread(cgname,'filtscale');
filtscale=filtscale(2:end);
Thm_Eulerian=Thm_Eulerian(2:end);

colors={'#0072BD','#D95319','#EDB120','#7E2F8E'};
colors_rgb = {...
    [0, 114/255, 189/255], ...   % #0072BD
    [217/255, 83/255, 25/255], ... % #D95319
    [237/255, 177/255, 32/255], ... % #EDB120
    [126/255, 47/255, 142/255] ... % #7E2F8E
};
kf1=fliplr(1./filtscale'.*2.*pi);
A={'a1','a2','a3','a4'}
B={'b1','b2','b3','b4'}

for ii=1:length(np)
    nparticles=np(ii);
    
    outputname=[Case,'_pars_P',num2str(nparticles),'T',num2str(days),...
        'bootstrap.mat']
    load(outputname)
    dstr=1;
    if strcmpi(Case, 'wave')
        if ii==1
            lambda=1e-9;
        else
            lambda=1e-10;
        end
    end

    if strcmpi(Case, 'nowave')
        if ii==1
            lambda=1e-8;
        elseif ii==2
            lambda=1e-9;
        else
            lambda=1e-10;
        end
    end

    % for tt=1:size(SF3,2)
    %      [SpecFlux(:,tt),Vt(:,tt),ebs(:,tt),kf,lf]=Fk_fitting_SF3(SF3(dstr:end,tt),dist_axis(dstr:end), ...
    %          2e3,500e3,'fuc','RLS',lambda,kf1);
    % end

    for tt=1:size(SF3,2)
         [SpecFlux(:,tt),Vt(:,tt),ebs(:,tt),kf,lf]=Fk_fitting_SF3(SF3(dstr:end,tt),dist_axis(dstr:end), ...
             2e3,500e3,'log','RLS',lambda,kf1);
    end
    
    % [residual_norms, solution_norms,...
    % lambda_opt_idx]=Fk_fitting_SF3_Lcurve(SF3_mean(dstr:end),dist_axis(dstr:end), ...
    % 1,300e3, 'fuc','RLS',[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
    %     1e-7,1e-8,1e-9,1e-10,1e-11,1e-12],kf1);
    clear CI_ebs CI_Vt CI_SpecFlux CI_SF3
    
    % Convert parameters to variance preserving form for plotting
    ebs_varp = ebs;
    ebs_varp(1:end-1,:) = ebs(1:end-1,:).*kf';
    
    for i = 1:size(ebs,1)
        CI_ebs(:,i) = prctile(ebs_varp(i,:), [95, 5]); 
    end
    
    for i = 1:size(Vt,1)
        CI_Vt(:,i) = prctile(Vt(i,:), [95,5]);
    end
    
    for i = 1:size(SF3,1)
        CI_SF3(:,i) = prctile(SF3(i,:), [95,5]);
    end
    
    for i = 1:size(SpecFlux,1)    
        CI_SpecFlux(:,i) = prctile(SpecFlux(i,:), [95,5]);
    end
    
    median_ebs = median(ebs_varp,2); 
    mean_ebs = nanmean(ebs_varp,2); 
    median_Vt = median(Vt,2);
    mean_Vt = nanmean(Vt,2);
    
    %median_S = median(S,2);
    mean_SF3 = nanmean(SF3,2);
    
    median_SpecFlux = median(SpecFlux,2);
    mean_SpecFlux = nanmean(SpecFlux,2);
    
    
    figure(1)
    semilogx(dist_axis./1e3,nanmean(SF3,2)./dist_axis','Color',colors{ii}, ...
    'LineWidth',1.5);
    hold on
    semilogx(lf./1e3,mean_Vt./lf','Color',colors{ii}, ...
    'LineWidth',1.5,'linestyle','--');
    grid on
    xlim([1e3,1e6]./1e3);
    xlabel('km')
    ylabel('m^{3}/s^{3}')
    title([Case,': sf3/r'])
    set(gca,'fontsize',16,'FontWeight','b')

    figure(2)
    A{ii}=semilogx(1./kf./1e3.*2.*pi,mean_SpecFlux,'Color',colors{ii}, ...
    'LineWidth',1.5);
    hold on
    x_fill = [1./kf./1e3.*2.*pi, fliplr(1./kf./1e3.*2.*pi)];
    y_fill = [CI_SpecFlux(1,:), fliplr(CI_SpecFlux(2,:))];
    % y_fill = [(nanmean(SpecFlux,2)+std1)', fliplr((nanmean(SpecFlux,2)-std1)')];
    fill(x_fill, y_fill, colors_rgb{ii}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

    eval(['load ',Case,'_pars_P',num2str(nparticles),'T89.5daysCG_Lag_spectukey', ...
    ini,'.mat']);
    
    Th_all=Th_all(2:end,:);
    Th=Th(2:end);
    std1=std(Th_all, 0, 2, 'omitnan');
    x_fill = [filtscale'./1e3, fliplr(filtscale'./1e3)];
    y_fill = [(Th'+std1)', fliplr((Th'-std1)')];
    B{ii}=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{ii},'LineStyle','--');
    fill(x_fill, y_fill, colors_rgb{ii}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    grid on
    xlim([1e3,1e6]./1e3);
    ylim([-0.4e-7,0.4e-7]);
    xlabel('km')
    ylabel('m^{2}/s^{3}')
    if ii==3
    eval(['load ',Case,'_Eulerian_SF3.mat'])
    N=287;dstr=1;dot=202;
    
    [r,SF3,S3L1,S3T1]=calc_radial(S3L1,S3T1,N,xscale);
    clear SpecFlux_E;clear kf_E;
    % [SpecFlux_E,Vt_E,ebs_E,kf_E,lf_E]=Fk_fitting_SF3(SF3(3:dot)', ...
    %     r(3:dot),1,300e3,'log','RLS',1e-10,kf1);
    [SpecFlux_E,Vt_E,ebs_E,kf_E,lf_E]=Fk_fitting_SF3(SF3(1:dot)', ...
        r(1:dot),1,300e3,'fuc','RLS',1e-10,kf1);
    c1=semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,(SpecFlux_E(dstr:end)),'LineWidth',1.5,'Color','k','LineStyle','-');
    % legend([A{1}, B{1}, A{2}, B{2}, A{3}, B{3},c1], ...
    % {'RLS-P289','CG-P289','RLS-P625','CG-P625', ...
    % 'RLS-P2500','CG-P2500','Eul fk'})
    end
    title([Case,': Moving Bootstrap RLS-fk'])
    set(gca,'fontsize',16,'FontWeight','b')

    % figure(3)
    % A{ii}=semilogx(1./kf./1e3,mean_ebs(1:end-1),'Color',colors{ii}, ...
    % 'LineWidth',1.5);
    % hold on
    % grid on

    figure(4)
    subplot(2,3,ii)
    A{ii}=semilogx(1./kf./1e3.*2.*pi,mean_SpecFlux,'Color',colors{ii}, ...
    'LineWidth',1.5);
    hold on
    
    B{ii}=semilogx(filtscale./1e3,Th', ...
        'LineWidth',1.5,'Color',colors{ii},'LineStyle','--');
    grid on
    ylim([-5e-8,5e-8]);
    

    figure(5)
    ebs1=nanmean(ebs,2);
    semilogx(kf.*1e3./(2.*pi),ebs1(1:end-1).*abs(mean(diff(kf))),'LineWidth',1.5, ...
    'Color',colors{ii});
    % abs(mean(diff(kf)))
    hold on
    grid on
end
