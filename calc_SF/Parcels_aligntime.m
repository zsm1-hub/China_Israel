clear
Case='wave'; % wave
nparticles=289*2148/6; % numbers of particles
days=89.5;  % days
seconds=0.2;
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
% input_dir='/meddy/simingzhang/Data/Parcels_data/';
% input_dir='/meddy/simingzhang/Data/Parcels_data/onetime_spectukey/';
input_dir='/meddy/simingzhang/Data/Parcels_data/';
addpath(genpath('/meddy/simingzhang/Analysis/matlab/Parcels_SF'))


% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
timerange=1:2140;

if strcmpi(Case, 'wave')
    fname=[input_dir,'wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
end

if strcmpi(Case, 'nowave')
    fname=[input_dir,'nowave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
end

if strcmpi(Case, 'HIT2d')
    fname=[input_dir,Case,'_pars_P',num2str(nparticles),'T',num2str(seconds),'seconds.nc'];
end

% clear
% fname='wave_pars_P103462T89.5days.nc'

info = ncinfo(fname);
var_names = {info.Variables.Name};
Time=ncread(fname,'time');
% nc=netcdf(fname,'w');

for ii=1:length(var_names)
    % disp()
    dt=3600;
    eval([var_names{ii},'=ncread(fname,','''',var_names{ii},'''',');']);
    eval(['[Taxis,Traxis]','=(size(',var_names{ii},'));'])
    eval([upper(var_names{ii}),'=zeros(',num2str(Taxis),',',num2str(Traxis),');'])
    eval([upper(var_names{ii}),'=trans_time(',var_names{ii},',dt,Time);'])
    % LAT=trans_time(lat,dt,Time);
    eval(['clear ',var_names{ii}]);
    eval([var_names{ii},'=',upper(var_names{ii}),';']);
    eval(['clear ',upper(var_names{ii})]);

    % eval(['nc{','''',var_names{ii},'''','}=',var_names{ii},'''',';']);

    disp([var_names{ii}, ' done']);
end

save([fname(1:end-3),'traj.mat'],'lat','lon','ue','ve');

% close(nc)