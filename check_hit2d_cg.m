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
a1=semilogx(1./K1D,specFlux_mean,'LineWidth',1.5, ...
'Color','k','LineStyle','-'); 
hold on

%2. Eul cg spec
fname='s2sflux_spec_hit_tukey.0002.nc';
ncdisp(fname)
Thm_Eulerian=ncread(fname,'Thm');
filtscale=ncread(fname,'filtscale');
a2=semilogx(filtscale,Thm_Eulerian,'LineWidth',1.5,'Color', ...
    'b');


%3. Eul cg uni
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
[SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(SF3(1:dot)',r(1:dot),0.009,6.32, ...
    'fuc','RLS',lambda,kf1);
% [residual_norms, solution_norms,...
% lambda_opt_idx]=Fk_fitting_SF3_Lcurve(SF3(1:dot)',r(1:dot),0.009,6.32, ...
%     'fuc','RLS',[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
%     1e-7,1e-8,1e-9,1e-10,1e-11,1e-12],kf1);
a4=semilogx(1./kf.*(2.*pi),SpecFlux,'LineWidth',1.5, ...
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
% legend([a1,a2,a3,a4],{'FFT specflux','CG (sfilt)','CG (unifilt)','RLS (angular k)'})
legend([a1,a2,a3,a4], ...
    {'FFT specflux','CG (sfilt)','CG (unifilt)','RLS (angular k)'}, ...
    "location",'northeast')
title(['HIT2d: CG vs RLS cross-scale energy flux'],'Interpreter','latex')
set(gca,'fontsize',14,'FontWeight','b')

% % subregion
% inset_pos = [0.18, 0.18, 0.18, 0.18]; % [left, bottom, width, height]
% ax_inset = axes('Position', inset_pos);
% box on; % 添加边框
% semilogx(1./K1D,specFlux_mean,'LineWidth',1.5, ...
% 'Color','magenta','LineStyle','-'); 
% hold on
% 
% fname='s2sflux_spec_hit_tukey.0002.nc';
% ncdisp(fname)
% Thm_Eulerian=ncread(fname,'Thm');
% filtscale=ncread(fname,'filtscale');
% semilogx(filtscale,Thm_Eulerian,'LineWidth',1.5,'Color', ...
%     colors{jj});
% 
% fname='s2sflux_spec_hit_uni.0002.nc';
% ncdisp(fname)
% Thm_Eulerian=ncread(fname,'Thm');
% filtscale=ncread(fname,'filtscale');
% semilogx(filtscale(2:end),Thm_Eulerian(2:end),'LineWidth',1.5,'Color','k');
% 
% load HIT2d_Eulerian_SF3.mat
% lambda=8e-1;
% dot=362;
% 
% %5. use kf1 to fit
% kf1=K1D.*2.*pi;
% % kf1=1./filtscale.*2.*pi;
% [SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(SF3(1:dot)',r(1:dot),0.009,6.32, ...
%     'fuc','RLS',lambda,kf1);
% semilogx(1./kf.*(2.*pi),SpecFlux,'LineWidth',1.5, ...
% 'Color','r');
% grid on
% xlim([2*pi/512*2,1e-1])
% ylim([-3,3])

%% energy injection
figure(8)
% dkk=mean(abs(K1D(2:end)-K1D(1:end-1)));
semilogx(K1D.*2.*pi,-divFlux_mean,'LineWidth',1.5, ...
'Color','k')
hold on
grid on
% dk=-(kf(2:end)-kf(1:end-1)).*2.*pi
semilogx(kf,ebs(1:end-1).*abs(mean(diff(kf))),'LineWidth',1.5, ...
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
dPidk_s=(Thmc(2:end)-Thmc(1:end-1))./(kc(2:end)-kc(1:end-1))
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
dPidk_u=(Thmc(2:end)-Thmc(1:end-1))./(kc(2:end)-kc(1:end-1))

semilogx(kc_mid_u.*2.*pi,dPidk_u,'LineWidth',1.5,'Color', ...
    [0.7, 0.9, 0.7]);

% plot([40.6,25].*2.*pi,[0,-37],'LineWidth',1.0,'Color','k')
% plot([4.2,0.77].*2.*pi,[0,-37],'LineWidth',1.0,'Color','k')
plot([40.6,40].*2.*pi,[0,-37],'LineWidth',1.0,'Color','k')
plot([4.2,1.2].*2.*pi,[0,-37],'LineWidth',1.0,'Color','k')

ylim([-120,10])
xlabel('k [1/m]')
ylabel('\epsilon_{j} dk_{j} [m^{2}/s^{3}]')
% legend([a1,a2,a3,a4,a5],{'cg spec','cg uni','fft-specflux','RLS','RLS (angular k)'})
% legend([a1,a2,a3,a4],{'FFT specflux','CG (sfilt)','CG (unifilt)','RLS (angular k)'})
legend([a1,a2,a3,a4], ...
    {'FFT specflux','CG (sfilt)','CG (unifilt)','RLS (angular k)'}, ...
    "location",'northeast')
title(['HIT2d: CG vs RLS energy injection'],'Interpreter','latex')
set(gca,'fontsize',14,'FontWeight','b')

% subregion
inset_pos = [0.38, 0.28, 0.38, 0.38]; % [left, bottom, width, height]
ax_inset = axes('Position', inset_pos);
box on; % 添加边框

semilogx(K1D.*2.*pi,-divFlux_mean,'LineWidth',1.5, ...
'Color','k');hold on
semilogx(kc_mid_s.*2.*pi,dPidk_s,'LineWidth',1.5,'Color', ...
    'b');
semilogx(kc_mid_u.*2.*pi,dPidk_u,'LineWidth',1.5,'Color', ...
    [0.7, 0.9, 0.7]);
semilogx(kf,ebs(1:end-1).*abs(mean(diff(kf))),'LineWidth',1.5, ...
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

%%

figure(2)
a1=semilogx(r,(SF3)./r,'Color',colors{2},'LineWidth',1.5)
% semilogx(dist_axis./1e3,SF1L,'Marker','x','Color',colors{ii},'LineWidth',1.5,'LineStyle','--')
% semilogx(dist_axis./1e3,SF1T,'Marker','x','Color',colors{ii},'LineWidth',1.5,'LineStyle','-.')
hold on
a2=semilogx(lf,Vt./lf','Color','r','LineWidth',1.5,'LineStyle','none', ...
    'marker','+')

% loglog(r,-(SF3),'Color',colors{ii},'LineWidth',1.5,'Marker','+')
% loglog(r,(SF3),'Color',colors{ii},'LineWidth',1.5,'Marker','o')
% [x32,y32]=get_line_loglog(3/2,1e-2,1e-2,-2,0);
% loglog(x32,y32,'LineWidth',1.5,'color','k')
% [x45,y45]=get_line_loglog(4/5,1e-2,1e-2,-2,0);
% loglog(x45,y45,'LineWidth',1.5,'color','b')
grid on
xlabel('r [m]')
ylabel('SF3/r [m^{2}/s^{3}]')
% legend([a1,a2,a3,a4,a5],{'cg spec','cg uni','fft-specflux','RLS','RLS (angular k)'})
% legend([a1,a2,a3,a4],{'FFT specflux','CG (sfilt)','CG (unifilt)','RLS (angular k)'})
legend([a1,a2], ...
    {'Third-order structure function','RLS-fitting'}, ...
    "location",'northeast')
title(['HIT2d: Third-order structure function'],'Interpreter','latex')
set(gca,'fontsize',14,'FontWeight','b')


% figure(4)
% semilogx(r,(SF3)./r,'Color',colors{ii},'LineWidth',1.5);hold on
% semilogx(r,nanmean(SF3_time,1)./r,'Color',colors{ii+1},'LineWidth',1.5, ...
%     'Marker','+','LineStyle','none');

% figure(5)
% for j=1:199
%     [SpecFlux(:,j),Vt(:,j),ebs(:,j),kf,lf]=Fk_fitting_SF3(SF3_time(j,1:dot)', ...
%         r(1:dot)./sqrt(2),0.009,6.32,'log','RLS',lambda);
% end
% std1=std(SpecFlux, 0, 2, 'omitnan');
% x_fill = [1./kf.*rescale, fliplr(1./kf.*rescale)];
% % y_fill = [CI_SpecFlux(1,:), fliplr(CI_SpecFlux(2,:))];
% y_fill = [(nanmean(SpecFlux,2)+std1)', fliplr((nanmean(SpecFlux,2)-std1)')];
% semilogx(1./kf.*rescale,nanmean(SpecFlux,2),'Color',colors{ii},'LineWidth',1.5);hold on
% fill(x_fill, y_fill, colors_rgb{ii}, 'FaceAlpha', 0.3, 'EdgeColor', 'none')
% 
% semilogx(filtscale,Thm_Eulerian,'LineWidth',1.5,'Color',colors{jj});hold on
% fname='s2sflux_spec_hit_uni.0002.nc';
% ncdisp(fname)
% Thm_Eulerian=ncread(fname,'Thm');
% filtscale=ncread(fname,'filtscale');
% semilogx(filtscale,Thm_Eulerian,'LineWidth',1.5,'Color','k');
% 
% fname='s2sflux_spec_hit_tukey.0002.nc';
% ncdisp(fname)
% Thm_Eulerian=ncread(fname,'Thm');
% filtscale=ncread(fname,'filtscale');
% semilogx(filtscale,Thm_Eulerian,'LineWidth',1.5,'Color',[.7,.7,.7]);
% grid on
% ylim([-20,20])
% 
% figure(6)
% semilogx(r,(S3T1_time+S3L1_time)','LineWidth',1.5);hold on
% grid on


% load HIT2d_pars_P2500T0.2secondsCG_Lag_uni_grid.mat
% a3=semilogx(filtscale,Th','LineWidth',1.5,'Color',colors{ii+1},'LineStyle','--');

% a7=semilogx(1./kf(dstr:end)./1e3.*rescale,SpecFlux(dstr:end),'LineWidth',1.5, ...
% 'Color','k','LineStyle','--');

% nparticles=[15376];
% clear fname
% fname{ii,:}=[Case,'_pars_P',num2str(nparticles(ii)),'T',num2str(89.5),'days.nc'];
% eval(['load ',fname{ii}(1:end-3),'SF123',ini,'.mat']);
% if ii>2
%     SF3=(SF3lll_time+SF3ltt_time);
%     SF1=(SF1l_time+SF1t_time)';
%     SF1L=(SF1l_time)';
%     SF1T=(SF1t_time)';
% else
%     SF3=(SF3lll+SF3ltt)';
%     SF1=(SF1l+SF1t)';
%     SF1L=(SF1l)';
%     SF1T=(SF1t)';
% 
% end
% [SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(SF3,dist_axis,1,300e3,'log','RLS',lambda);
% dstr=14;
% % dstr=1;
% a4=semilogx(1./kf(dstr:end)./1e3,SpecFlux(dstr:end),'Marker','x','Color',colors{ii+1},'LineWidth',1.5);
% eval(['load ',Case,'_pars_P',num2str(nparticles(ii)),'T89.5daysCG_Lag_spectukey',ini,'.mat']);
% a2=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{ii+1},'LineStyle','--');
% 
% grid on
% legend(['a1','a2','a4','a6'],{'Eulerian cg fk','Eulerian-SF3-fitting fk', ...
%     'P15376 Lag SF3-fitting fk','P15376 Lag cg fk'}, ...
%     'Location', 'northeast')
% ylim([-0.2e-7,0.5e-7]);
% xlim([1e3,1e6]./1e3);
% xlabel('km')
% ylabel('m^{2}/s^{3}')
% title('Smooth case')
% set(gca,'fontsize',16,'FontWeight','b')
