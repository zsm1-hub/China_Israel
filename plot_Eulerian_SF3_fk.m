
%%
clear all;close all
clc
%%%%%%%%%%%%%%%%calc Eulerian SF3 for Iceland%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/meddy/simingzhang/Analysis/matlab/Parcels_SF/')
addpath('/meddy/simingzhang/Data/Parcels_data')
addpath('/meddy/simingzhang/Data/RB_iceland_data')
input_dir='/meddy/simingzhang/Data/RB_iceland_data/';

colors={'#0072BD','#D95319','#EDB120','#7E2F8E'};

%%%%%%% wave
Case='wave';
wavename=[Case,'_Eulerian_SF3.mat'];
eval(['load ',wavename])
N=287;
[r,SF3,S3L1,S3T1]=calc_radial(S3L1,S3T1,N,xscale);


figure(1)
ii=1
loglog(r,abs(SF3),'LineWidth',1.5,'Color',colors{ii});hold on
loglog(r,-(SF3),'LineWidth',1.5,'Color',colors{ii},'Marker','+')
loglog(r,(SF3),'LineWidth',1.5,'Color',colors{ii},'Marker','o')

dot=202
lambda=1e-10;
[SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(SF3(1:dot)',r(1:dot),1,200e3,'log','RLS',lambda);

figure(2)
ii=1
semilogx(1./kf./1e3,SpecFlux,'LineWidth',1.5,'Color',colors{ii});hold on
%%%%%%%%%nowave
Case='nowave';
wavename=[Case,'_Eulerian_SF3.mat'];
eval(['load ',wavename])
N=287;
[r,SF3,S3L1,S3T1]=calc_radial(S3L1,S3T1,N,xscale);


ii=2
figure(1)
loglog(r,abs(SF3),'LineWidth',1.5,'Color',colors{ii});hold on
loglog(r,-(SF3),'LineWidth',1.5,'Color',colors{ii},'Marker','+')
loglog(r,(SF3),'LineWidth',1.5,'Color',colors{ii},'Marker','o')


dot=202
lambda=1e-10;
[SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(SF3(1:dot)',r(1:dot),1,200e3,'log','RLS',lambda);

figure(2)
ii=2
semilogx(1./kf./1e3,SpecFlux,'LineWidth',1.5,'Color',colors{ii});hold on
grid on


