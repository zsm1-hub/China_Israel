%% check manuel code
clear all;close all;clc
fname='mSF_15.nc'
% ncdisp(fname)
load fftspecflux.mat

% semilogx(K1D,specFlux_mean)
SF3=ncread(fname,'du3');
r=nanmean(ncread(fname,'dr'),2);
semilogx(r,nanmean(SF3,2))

lambda=1e-10;
[SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(nanmean(SF3,2), ...
    r',r(1),r(end-1),'fuc','RLS',lambda,K1D.*2.*pi);

figure(2)
semilogx(K1D,specFlux_mean)
hold on
semilogx(kf./(2.*pi),SpecFlux)

