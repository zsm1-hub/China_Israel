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
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
input_dir='/meddy/simingzhang/Data/Parcels_data/';
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
% timerange=1:2140;
inv_style='RLS';
lambda=1e-10;

if strcmpi(Case, 'wave')
    fname=['wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
end

if strcmpi(Case, 'nowave')
    fname=['nowave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
end

eval(['load ',Case,'_pars_P',num2str(nparticles),'T',num2str(days),'daysSF123.mat']);

%% Fitting NNLS

SF3=(SF3lll_time+SF3ltt_time)';
nsamps = size(SF3,2);

r=dist_axis;
Nr=length(r);

% select part of the r axis that we think has reasonable
% data.
% modified by zsm
% ns=find(dist_axis>=2e3 ,1);
% ne=find(dist_axis<=500e3 ,1,'last');

ns=find(dist_axis>=2 ,1);
ne=find(dist_axis<=500e3 ,1,'last');

R=r(ns:ne); 

NR=length(R);
lf=R;


kf = logspace( log10(1/max(R)), log10(1/min(R)), length(R)-1) ;
dk = diff(kf);
kf = 0.5*(kf(1:end-1) + kf(2:end));

% Code is written in a way that it is assumed that
% kf is decreasing. This matters when doing cumsum
kf = fliplr(kf); 
dk = fliplr(dk);

Nk=length(kf);
%%

% Define matrices
ebs = zeros(Nk+1,nsamps);
S = zeros(Nr, nsamps);
Vt = zeros(NR,nsamps);
SpecFlux = zeros(Nk, nsamps); 

norm_flag=1;

for n=1:nsamps
    
    % SF3 from the data
    %S(:,n) =s3lll(:,n)' +s3ltt(:,n)';
    S(:,n) = SF3(:,n);
    
    % SF3 over the selected range od points
    V=S(ns:ne, n)';
    
    % We will solve Ax=b problem
    % where x is the parameters (energy fluxes), b is SF3
    
    % build the matrix A
    A=zeros(NR,Nk+1);
    for j=1:Nk
        A(:,j)=-4/kf(j)*besselj(1,kf(j)*R)'*dk(j);
    end
    A(:,end)=2*R';
    
    % %% divide Ax=b by r on both sides to remove log dependence
    % normalize the magnitude
    if norm_flag == 1
        for j=1:NR
            A(j,:)=A(j,:)./abs(R(j));
        end
        V=V./abs(R);
    end
    
    % estimate the epsilons using least squares
    if strcmp(inv_style, 'NNLS')
        ebs(:,n) = lsqnonneg(A,V');
    elseif strcmp(inv_style, 'LS')
        ebs(:,n) = A\(V');

    elseif strcmp(inv_style, 'RLS')
        % 确保V是列向量
        V = V(:);  % 将行向量转为列向量
        
        % 获取A的列数（未知数个数）
        n_cols = size(A, 2);  % 应该是Nk+1 (27)
        
        % 构建正则化系统
        A_aug = [A; sqrt(lambda) * eye(n_cols)];
        V_aug = [V; zeros(n_cols, 1)];
        
        % 求解并确保结果是列向量
        ebs(:, n) = A_aug \ V_aug;
    end
    
    % What does the reconstructed V look like
    Vt(:,n) =2*ebs(end,n)*R;
    for j=1:Nk
        Vt(:,n)=Vt(:,n)-4*ebs(j,n)./kf(j).*besselj(1,kf(j).*R)'*dk(j);
    end
    
    % Convert V to spectral fluxes   
    SpecFlux(1,n) = -ebs(end,n); 
    
    for j = 2:Nk
        SpecFlux(j,n) = SpecFlux(j-1,n) + ebs(end-j+1, n)*dk(end-j+2);
    end
    
end



subplot(222)
semilogx(dist_axis./1e3,SF3'); hold on
semilogx(lf./1e3,Vt)
title('SF3')
xlabel('km')
ylabel('m^{2}/s^{3}')
grid on


subplot(221)
SpecFlux = flipud(SpecFlux); 
semilogx(2*pi.*1./kf./1e3,SpecFlux,'Marker','x')
title('F(k)')
xlabel('km')
ylabel('m^{2}/s^{3}')
grid on



subplot(223)
semilogx(2*pi./kf./1e3,ebs(1:end-1),'Marker','x')
title('ebs')
xlabel('km')
grid on

