%% check hit2d
clear;close all;clc

% lambda=1e-10;
rescale=1
ii=4
jj=1
colors={[.7,.7,.7],'#0072BD','#D95319','#EDB120','#7E2F8E'};
colors_rgb = {...
    [.7,.7,.7],...
    [0, 114/255, 189/255], ...   % #0072BD
    [217/255, 83/255, 25/255], ... % #D95319
    [237/255, 177/255, 32/255], ... % #EDB120
    [126/255, 47/255, 142/255] ... % #7E2F8E
};

%% CG vs RLS cross-scale energy flux
figure(1)

%1. fftspec
load HIT2d_fftfk.mat
a1=semilogx(1./K1D,specFlux_mean,'LineWidth',2, ...
'Color','k','LineStyle','-'); 
hold on

%2. Eul cg spec
fname='s2sflux_spec_hit_tukey.0002.nc';
ncdisp(fname)
Thm_Eulerian=ncread(fname,'Thm');
filtscale=ncread(fname,'filtscale');
a2=semilogx(filtscale,Thm_Eulerian,'LineWidth',2,'Color', ...
    'b');


% 3. Eul cg uni
fname='s2sflux_spec_hit_uni.0002.nc';
ncdisp(fname)
Thm_Eulerian=ncread(fname,'Thm');
filtscale=ncread(fname,'filtscale');
a3=semilogx(filtscale(2:end),Thm_Eulerian(2:end),'LineWidth',1.5,'Color', ...
    [0.7, 0.9, 0.7]);


%4. use log fit and didn't plus 2pi
load HIT2d_Eulerian_SF3.mat
lambda=8e-1;
dot=362;

% kf1=fliplr(1./filtscale(2:end)'.*2.*pi);
% [SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(SF3(1:dot)',r(1:dot),0.009,6.32, ...
%     'log','RLS',lambda,kf1);
% 
% a4=semilogx(1./kf.*rescale,SpecFlux,'LineWidth',1.5, ...
% 'Color',colors{jj+1});

%5. use kf1 to fit
kf1=K1D.*2.*pi;
% kf1=1./filtscale.*2.*pi;
[SpecFlux_E,Vt_E,ebs_E,kf_E,lf_E]=Fk_fitting_SF3(SF3(1:dot)',r(1:dot),0.009,6.32, ...
    'fuc','RLS',lambda,kf1);


a4=semilogx(1./kf_E.*(2.*pi),SpecFlux_E,'LineWidth',2, ...
'Color','r');

text(1/64*2*pi+0.02,-0.9,'Injection scale: $\frac{1}{64} \times 2\pi \approx 0.1$', 'Interpreter', 'latex', ...
    'FontSize',14,'FontWeight','b')
plot([1/64*2*pi,1/64*2*pi],[-20,10],'LineStyle','--','color','k', ...
    'LineWidth',1.2)

grid on
ylim([-10,5]);
xlabel('r [m]')
ylabel('\Pi [m^{2}/s^{3}]')

% legend([a1,a2,a3,a4,a5],{'cg spec','cg uni','fft-specflux','RLS','RLS (angular k)'})
legend([a1,a2,a3,a4],{'FFT specflux','CG (sfilt)','CG (unifilt)','RLS (angular k)'})
% legend([a1,a2,a3,a4,a5], ...
%     {'FFT specflux','CG (sfilt)','CG (unifilt)','RLS (angular k)',['RLS ', ...
%     num2str(timerange(1)),'~',num2str(timerange(end))]}, ...
%     "location",'northeast')
title(['Eul: CG vs RLS cross-scale energy flux'],'Interpreter','latex')
set(gca,'fontsize',14,'FontWeight','b')

%%
figure(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%load Lag SF3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str1=0.05;
en1=5;
% lambda_opt_idx=[100,10,1,1e-1,1e-2];
lambda_opt_idx=[8e-1,8e-2,8e-3,8e-4];

load HIT2d_pars_P10000T150timeavg4.mat
r_10000=dist_axis;
SF3_10000=SF3_mean;
Th_10000=nanmean(Th_all,2);
filtscale_L=filtscale(2:end);

[SpecFlux_L10000,Vt_L10000,ebs_L10000,...
    kf_L,lf_L]=Fk_fitting_SF3_Lcurve(SF3_10000, ...
    r_10000,str1,en1, ...
    'log','RLS',lambda_opt_idx,kf1,0);

load HIT2d_pars_P62500T150timeavg4.mat
r_62500=dist_axis;
SF3_62500=SF3_mean;
Th_62500=nanmean(Th_all,2);

[SpecFlux_L62500,Vt_L62500,ebs_L62500,...
    kf_L,lf_L]=Fk_fitting_SF3_Lcurve(SF3_62500, ...
    r_62500,str1,en1, ...
    'log','RLS',lambda_opt_idx,kf1,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a1=semilogx(filtscale_L,Th_10000,'LineWidth',1.5, ...
'Color', [0.8, 0.5, 0.2],'LineStyle','--');
hold on
a2=semilogx(filtscale_L,Th_62500,'LineWidth',1.5, ...
'Color',[0.5, 0.6, 0.4] ,'LineStyle','--');

a3=semilogx(1./kf_L.*(2.*pi),(SpecFlux_L10000),'LineWidth',1.5, ...
'Color', [0.8, 0.5, 0.2],'LineStyle','-');

a4=semilogx(1./kf_L.*(2.*pi),(SpecFlux_L62500),'LineWidth',1.5, ...
'Color',[0.5, 0.6, 0.4] ,'LineStyle','-');
plot([1/64*2*pi,1/64*2*pi],[-20,10],'LineStyle','--','color','k', ...
    'LineWidth',1.2)

grid on
xlim([1e-2,1e1]);
ylim([-10,5]);
xlabel('r [m]')
ylabel('\Pi [m^{2}/s^{3}]')


legend([a1,a2,a3,a4],{'Lag CG P10000','Lag CG P62500','Lag RLS P10000','Lag RLS P62500'})
title(['Lag: CG vs RLS cross-scale energy flux'],'Interpreter','latex')
set(gca,'fontsize',14,'FontWeight','b')
%% Eul energy injection
figure(3)
% dkk=mean(abs(K1D(2:end)-K1D(1:end-1)));
semilogx(K1D.*2.*pi,-divFlux_mean,'LineWidth',1.5, ...
'Color','k')
hold on
grid on
% dk=-(kf(2:end)-kf(1:end-1)).*2.*pi
semilogx(kf_E,ebs_E(1:end-1).*abs(mean(diff(kf_E))),'LineWidth',1.5, ...
'Color','r');
% semilogx(kf./(2.*pi),ebs(1:end-1).*kf'./(2.*pi),'LineWidth',1.5, ...
% 'Color','r');

% for cg dPi/dk=-T(k)
%2. Eul cg spec
fname='s2sflux_spec_hit_tukey.0002.nc';
ncdisp(fname)
Thm_Eulerian=ncread(fname,'Thm');
filtscale=ncread(fname,'filtscale');
kc=flipud(1./filtscale);
Thmc=flipud(Thm_Eulerian);
kc_mid_s=0.5.*(kc(2:end)+kc(1:end-1));
dPidk_s=(Thmc(2:end)-Thmc(1:end-1))
% ./(kc(2:end)-kc(1:end-1))
% a2=semilogx(filtscale,Thm_Eulerian,'LineWidth',1.5,'Color', ...
%     'b');
semilogx(kc_mid_s.*2.*pi,dPidk_s,'LineWidth',1.5,'Color', ...
    'b');

fname='s2sflux_spec_hit_uni.0002.nc';
ncdisp(fname)
Thm_Eulerian=ncread(fname,'Thm');
filtscale=ncread(fname,'filtscale');
filtscale=filtscale(2:end);
Thm_Eulerian=Thm_Eulerian(2:end);
kc=flipud(1./filtscale);
Thmc=flipud(Thm_Eulerian);
kc_mid_u=0.5.*(kc(2:end)+kc(1:end-1));
dPidk_u=(Thmc(2:end)-Thmc(1:end-1))
% ./(kc(2:end)-kc(1:end-1))
semilogx(kc_mid_u.*2.*pi,dPidk_u,'LineWidth',1.5,'Color', ...
    [0.7, 0.9, 0.7]);

% plot([40.6,25].*2.*pi,[0,-37],'LineWidth',1.0,'Color','k')
% plot([4.2,0.77].*2.*pi,[0,-37],'LineWidth',1.0,'Color','k')
plot([40.1,40].*2.*pi,[0,-4.7],'LineWidth',1.0,'Color','k')
plot([4.2,1.2].*2.*pi,[0,-4.7],'LineWidth',1.0,'Color','k')

% ylim([-120,10])
ylim([-20,10])

xlabel('k [1/m]')
ylabel('\epsilon_{j} dk_{j} [m^{2}/s^{3}]')
% legend([a1,a2,a3,a4,a5],{'cg spec','cg uni','fft-specflux','RLS','RLS (angular k)'})
% legend([a1,a2,a3,a4],{'FFT specflux','CG (sfilt)','CG (unifilt)','RLS (angular k)'})
% legend([a1,a2,a3,a4], ...
%     {'FFT specflux','CG (sfilt)','CG (unifilt)','RLS (angular k)'}, ...
%     "location",'northeast')
title(['Eul: CG vs RLS energy injection'],'Interpreter','latex')
set(gca,'fontsize',14,'FontWeight','b')

% subregion
inset_pos = [0.37, 0.18, 0.38, 0.38]; % [left, bottom, width, height]
ax_inset = axes('Position', inset_pos);
box on; % 添加边框

semilogx(K1D.*2.*pi,-divFlux_mean,'LineWidth',1.5, ...
'Color','k');hold on
semilogx(kc_mid_s.*2.*pi,dPidk_s,'LineWidth',1.5,'Color', ...
    'b');
semilogx(kc_mid_u.*2.*pi,dPidk_u,'LineWidth',1.5,'Color', ...
    [0.7, 0.9, 0.7]);
semilogx(kf_E,ebs_E(1:end-1).*abs(mean(diff(kf_E))),'LineWidth',1.5, ...
'Color','r');
grid on
plot([64,64],[-2,8],'LineStyle','--','color','k', ...
    'LineWidth',1.2)
text(64-30,6,'$k=64$', 'Interpreter', 'latex', ...
    'FontSize',14,'FontWeight','b')
ylim([-2,8])
xlim([4.2,40.6].*2.*pi)
xticklabels('')
set(gca,'fontsize',14,'FontWeight','b')

%% Lag energy injection
figure(4)
% dkk=mean(abs(K1D(2:end)-K1D(1:end-1)));
semilogx(K1D.*2.*pi,-divFlux_mean,'LineWidth',1.5, ...
'Color','k')
hold on
grid on
% dk=-(kf(2:end)-kf(1:end-1)).*2.*pi
semilogx(kf_L,ebs_L10000(1:end-1).*abs(mean(diff(kf_L))),'LineWidth',1.5, ...
'Color',[0.8, 0.5, 0.2]);

semilogx(kf_L,ebs_L62500(1:end-1).*abs(mean(diff(kf_L))),'LineWidth',1.5, ...
'Color',[0.5, 0.6, 0.4]);

% plot([40.6,25].*2.*pi,[0,-37],'LineWidth',1.0,'Color','k')
% plot([4.2,0.77].*2.*pi,[0,-37],'LineWidth',1.0,'Color','k')
plot([40.1,40].*2.*pi,[0,-4.7],'LineWidth',1.0,'Color','k')
plot([4.2,1.2].*2.*pi,[0,-4.7],'LineWidth',1.0,'Color','k')

% ylim([-120,10])
ylim([-20,10])

xlabel('k [1/m]')
ylabel('\epsilon_{j} dk_{j} [m^{2}/s^{3}]')
% legend([a1,a2,a3,a4,a5],{'cg spec','cg uni','fft-specflux','RLS','RLS (angular k)'})
% legend([a1,a2,a3,a4],{'FFT specflux','CG (sfilt)','CG (unifilt)','RLS (angular k)'})
% legend([a1,a2,a3,a4], ...
%     {'FFT specflux','CG (sfilt)','CG (unifilt)','RLS (angular k)'}, ...
%     "location",'northeast')
title(['Lag: CG vs RLS energy injection'],'Interpreter','latex')
set(gca,'fontsize',14,'FontWeight','b')

% subregion
inset_pos = [0.37, 0.18, 0.38, 0.38]; % [left, bottom, width, height]
ax_inset = axes('Position', inset_pos);
box on; % 添加边框

semilogx(K1D.*2.*pi,-divFlux_mean,'LineWidth',1.5, ...
'Color','k');hold on
semilogx(kf_L,ebs_L10000(1:end-1).*abs(mean(diff(kf_L))),'LineWidth',1.5, ...
'Color',[0.8, 0.5, 0.2]);

semilogx(kf_L,ebs_L62500(1:end-1).*abs(mean(diff(kf_L))),'LineWidth',1.5, ...
'Color',[0.5, 0.6, 0.4]);

grid on
plot([64,64],[-2,8],'LineStyle','--','color','k', ...
    'LineWidth',1.2)
text(64-30,6,'$k=64$', 'Interpreter', 'latex', ...
    'FontSize',14,'FontWeight','b')
ylim([-2,8])
xlim([4.2,40.6].*2.*pi)
xticklabels('')
set(gca,'fontsize',14,'FontWeight','b')

%% Eul SF3/3 and fit

figure(5)
a1=semilogx(r,(SF3)./r,'Color',colors{2},'LineWidth',1.5)
% semilogx(dist_axis./1e3,SF1L,'Marker','x','Color',colors{ii},'LineWidth',1.5,'LineStyle','--')
% semilogx(dist_axis./1e3,SF1T,'Marker','x','Color',colors{ii},'LineWidth',1.5,'LineStyle','-.')
hold on
a2=semilogx(lf_E,Vt_E./lf_E','Color',colors{2},'LineWidth',1.5,'LineStyle','none', ...
    'marker','+')

% semilogx(r_10000,SF3_10000./r_10000','LineWidth',1.5,'Color',[0.8, 0.5, 0.2])
% semilogx(lf_L,Vt_L10000./lf_L','Linestyle', ...
%     'none','Color',[0.8, 0.5, 0.2],'Marker','+','LineWidth',1.5)
% 
% semilogx(r_62500,SF3_62500./r_62500','LineWidth',1.5,'Color',[0.5, 0.6, 0.4])
% semilogx(lf_L,Vt_L62500./lf_L','Linestyle', ...
%     'none','Color',[0.5, 0.6, 0.4],'Marker','+','LineWidth',1.5)
% 
% semilogx(kf_L,ebs_L10000(1:end-1).*abs(mean(diff(kf_L))),'LineWidth',1.5, ...
% 'Color',[0.8, 0.5, 0.2]);
% semilogx(kf_L,ebs_L62500(1:end-1).*abs(mean(diff(kf_L))),'LineWidth',1.5, ...
% 'Color',[0.5, 0.6, 0.4]);

% loglog(r,-(SF3),'Color',colors{ii},'LineWidth',1.5,'Marker','+')
% loglog(r,(SF3),'Color',colors{ii},'LineWidth',1.5,'Marker','o')
% [x32,y32]=get_line_loglog(3/2,1e-2,1e-2,-2,0);
% loglog(x32,y32,'LineWidth',1.5,'color','k')
% [x45,y45]=get_line_loglog(4/5,1e-2,1e-2,-2,0);
% loglog(x45,y45,'LineWidth',1.5,'color','b')
grid on
xlabel('r [m]')
ylabel('$$\mathbf{D3(r)/r \ [m^{2}/s^{3}]}$$','Interpreter','latex')
ylim([-10,30])
xlim([1e-2,1e1])
legend([a1,a2], ...
    {'$$\mathbf{D3(r)/r}$$','RLS-fitting'}, ...
    "location",'northeast','Interpreter','latex')
title(['Eul: Third-order structure function'],'Interpreter','latex')
set(gca,'fontsize',14,'FontWeight','b')


%% Lag SF3/3 and fit

figure(6)


a1=semilogx(r_10000,SF3_10000./r_10000','LineWidth',1.5,'Color',[0.8, 0.5, 0.2])
hold on

a2=semilogx(lf_L,Vt_L10000./lf_L','Linestyle', ...
    'none','Color',[0.8, 0.5, 0.2],'Marker','+','LineWidth',1.5)

a3=semilogx(r_62500,SF3_62500./r_62500','LineWidth',1.5,'Color',[0.5, 0.6, 0.4])
a4=semilogx(lf_L,Vt_L62500./lf_L','Linestyle', ...
    'none','Color',[0.5, 0.6, 0.4],'Marker','+','LineWidth',1.5)
a5=semilogx(r,(SF3)./r,'Color',colors{2},'LineWidth',1.5)


grid on
xlabel('r [m]')
ylim([-10,30])
xlim([1e-2,1e1])

ylabel('$$\mathbf{D3(r)/r \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% legend([a1,a2], ...
%     {'Third-order structure function','RLS-fitting'}, ...
%     "location",'northeast')
legend([a1,a2,a3,a4,a5], ...
    {'$$\mathbf{D3(r)/r \ P10000}$$','RLS-fitting',...
    '$$\mathbf{D3(r)/r \ P30000}$$','RLS-fitting','Eul D3(r)/r'}, ...
    "location",'northeast','Interpreter','latex')
title(['Lag: Third-order structure function'],'Interpreter','latex')
set(gca,'fontsize',14,'FontWeight','b')