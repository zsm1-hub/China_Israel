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
win='kaiser'
% ini='_rough_500m'
ini='_roughsmall_rot'

% str1=15;
% end1=28;
% str1=1;
% end1=19;
inv='RLS'
yy1=-0.6e-8;
yy2=-1.0e-8;
yy3=-1.4e-8;
yy4=-1.8e-8;

timerange=1:1940;
% param_Eul_cg_sf3fk_500m
% SF3_E=SF3_mean';
% r_E=r;

% 2km param
param_Eul_cg_sf3fk
SF3_E=SF3(2:end)';
r_E=r(2:end)';
ymin_fk=-3e-8;ymax_fk=6e-8;
ymin_sf=-1e-7;ymax_sf=1e-7;
% ymin_ebsdk=-3e-8;ymax_ebsdk=3e-8;
ymin_ebsdk=-2e-4;ymax_ebsdk=3e-4;
ymin_right = -4e-8;    % 
ymax_right = 6e-8;  % 
xmin=1;
xmax=1e3;

addpath('E:\DATA\RoyBarkan\RB_iceland_data\RB_iceland_data')
grdname='E:\DATA\RoyBarkan\RB_iceland_data\RB_iceland_data\niskin2km_500m_grd.nc';
lon_g=ncread(grdname,'lon_rho');
lat_g=ncread(grdname,'lat_rho');


load Eulerian_SF2_spec_Iceland.mat


if strcmpi(Case, 'wave')
    SF2_E=nanmean(SF2w,2);
    SF2ll_E=nanmean(S2lw,2);
    SF2tt_E=nanmean(S2tw,2);
    Casemean='Hf';
    if strcmpi(ini,'_rough')
        lambda=[1e-7,1e-8,1e-9];
        % str1=15;
        % end1=28;
        str1=1;
        end1=19;
        range1=2.5e3;
        range2=500e3;
    end
    if strcmpi(ini,'_roughsmall')
        % lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
        % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,1e-14];
        lambda=[1e-7,1e-8,1e-9,1e-10,1e-11];
        str1=1;
        end1=19;
        range1=2.5e3;
        range2=300e3;
    end
    if strcmpi(ini,'_roughsmall_500m')
         lambda=[1e-8];
        str1=1;
        end1=24;
        range1=1e3;
        range2=300e3;
    end
    if strcmpi(ini,'_rough_500m') || strcmpi(ini,'_roughbox100g_500m')
        lambda=[1e-7];
        str1=1;
        end1=24;
        % end1=17;
        range1=1e3;
        range2=300e3;
    end
   
    % range2=300e3;
end
if strcmpi(Case, 'nowave')
    SF2_E=nanmean(SF2nw,2);
    SF2ll_E=nanmean(S2lnw,2);
    SF2tt_E=nanmean(S2tnw,2);
    Casemean='Sm';
    % lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
    % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,1e-14,1e-15,1e-16];
    if strcmpi(ini,'_rough')
        % lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
        % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12];
        lambda=[1e-7,1e-8,1e-9];
        str1=1;
        end1=19;
        range1=2.5e3;
        range2=300e3;
    end
    if strcmpi(ini,'_roughsmall')
     %     lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
     % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,7e-14];
         lambda=[1e-7,1e-8,1e-9,1e-10,1e-11];
        str1=1;
        end1=18;
        range1=2.5e3;
        range2=300e3;
    end

    if strcmpi(ini,'_roughsmall_500m')
        lambda=[1e-8];
        str1=1;
        end1=24;
        range1=1e3;
        range2=300e3;
    end
    if strcmpi(ini,'_rough_500m') || strcmpi(ini,'_roughbox100g_500m')
         lambda=[5e-9];
        str1=1;
        % end1=24;
        end1=15;
        range1=1e3;
        range2=300e3;
    end
   

end



timerange=1:1940;

fname=[Case,'_pars_P',num2str(289),'T',num2str(timerange(end)),ini,'bootstrap.mat'];
load(fname);
SF2_L289=SF2;
SF2ll_L289=SF2ll;
SF2tt_L289=SF2tt;
tname=[Case,'_pars_P',num2str(289),'T',num2str(89.5),'days',ini,'traj.mat']
load(tname)
lon289=traj.trajmat_X;
lat289=traj.trajmat_Y;


fname=[Case,'_pars_P',num2str(625),'T',num2str(timerange(end)),ini,'bootstrap.mat'];
load(fname);
SF2_L625=SF2;
SF2ll_L625=SF2ll;
SF2tt_L625=SF2tt;
tname=[Case,'_pars_P',num2str(625),'T',num2str(89.5),'days',ini,'traj.mat']
load(tname)
lon625=traj.trajmat_X;
lat625=traj.trajmat_Y;


fname=[Case,'_pars_P',num2str(1089),'T',num2str(timerange(end)),ini,'bootstrap.mat'];
load(fname);
SF2_L1089=SF2;
SF2ll_L1089=SF2ll;
SF2tt_L1089=SF2tt;
tname=[Case,'_pars_P',num2str(1089),'T',num2str(89.5),'days',ini,'traj.mat']
load(tname)
lon1089=traj.trajmat_X;
lat1089=traj.trajmat_Y;

figure(1)
a1=loglog(rbin./1e3,SF2_E,'LineWidth',1.5,'Color','k');hold on
a2=loglog(dist_axis./1e3,SF2_L289/2,'LineWidth',1.5,'Color',colors_rgb{4});
a3=loglog(dist_axis./1e3,SF2_L625/2,'LineWidth',1.5,'Color',colors_rgb{5});
a4=loglog(dist_axis./1e3,SF2_L1089/2,'LineWidth',1.5,'Color',colors_rgb{6});
grid on
[x23,y23]=get_line_loglog(2/3,10^0.1,10^-2.5,0.1,0.8)
[x2,y2]=get_line_loglog(2,10^0.1,10^-2.5,0.1,0.8)
[x1,y1]=get_line_loglog(1,10^0.1,10^-2.5,0.1,0.8)

loglog(x23,y23,'LineWidth',1,'Color','k','LineStyle','--');
loglog(x1,y1,'LineWidth',1,'Color','k','LineStyle','--');
loglog(x2,y2,'LineWidth',1,'Color','k','LineStyle','--');
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{D2(r) \ [m^{2}/s^{2}]}$$','Interpreter','latex')
legend([a1,a2,a3,a4],{['Eul SF2'],['Lag P289'], ...
    ['Lag P625'],['Lag P1089']},'Location','southeast')
title('Eul vs Lag: SF2')
set(gca,'fontsize',14,'fontweight','b')

% loglog(rbin,SF2ll_E);hold on
% loglog(dist_axis,SF2ll_L289);
% loglog(dist_axis,SF2ll_L625);
% loglog(dist_axis,SF2ll_L1089);
% 
% loglog(rbin,SF2tt_E);hold on
% loglog(dist_axis,SF2tt_L289);
% loglog(dist_axis,SF2tt_L625);
% loglog(dist_axis,SF2tt_L1089);

figure(2)
hold on
plot(lon_g(1,:),lat_g(1,:),'LineWidth',1.5,'Color','r','LineStyle','--')
plot(lon_g(end,:),lat_g(end,:),'LineWidth',1.5,'Color','r','LineStyle','--')
plot(lon_g(:,1),lat_g(:,1),'LineWidth',1.5,'Color','r','LineStyle','--')
plot(lon_g(:,end),lat_g(:,end),'LineWidth',1.5,'Color','r','LineStyle','--')
% a1=scatter(lon1089(1,:),lat1089(1,:),10,'MarkerFaceColor', ...
%     colors_rgb{3},'MarkerEdgeColor','flat')
% a2=scatter(lon625(1,:),lat625(1,:),10,'MarkerFaceColor', ...
%     colors_rgb{2},'MarkerEdgeColor','flat')
a3=scatter(lon289(1,:),lat289(1,:),10,'MarkerFaceColor', ...
    colors_rgb{1},'MarkerEdgeColor','flat')
legend([a1,a2,a3],{'1089','625','289'})
xlabel('Longitude')
ylabel('Latitude')
title([Casemean,' initial domain'])
set(gca,'fontsize',14,'fontweight','b')


figure(3)
hold on
plot(lon_g(1,:),lat_g(1,:),'LineWidth',1.5,'Color','r','LineStyle','--')
plot(lon_g(end,:),lat_g(end,:),'LineWidth',1.5,'Color','r','LineStyle','--')
plot(lon_g(:,1),lat_g(:,1),'LineWidth',1.5,'Color','r','LineStyle','--')
plot(lon_g(:,end),lat_g(:,end),'LineWidth',1.5,'Color','r','LineStyle','--')
plot(lon289,lat289,'LineWidth',1.5)
xlabel('Longitude')
ylabel('Latitude')
title([Casemean,' P289'])
set(gca,'fontsize',14,'fontweight','b')

figure(4)
hold on
plot(lon_g(1,:),lat_g(1,:),'LineWidth',1.5,'Color','r','LineStyle','--')
plot(lon_g(end,:),lat_g(end,:),'LineWidth',1.5,'Color','r','LineStyle','--')
plot(lon_g(:,1),lat_g(:,1),'LineWidth',1.5,'Color','r','LineStyle','--')
plot(lon_g(:,end),lat_g(:,end),'LineWidth',1.5,'Color','r','LineStyle','--')
plot(lon625,lat625,'LineWidth',1.5)
xlabel('Longitude')
ylabel('Latitude')
title([Casemean,' P625'])
set(gca,'fontsize',14,'fontweight','b')

figure(5)
hold on
plot(lon_g(1,:),lat_g(1,:),'LineWidth',1.5,'Color','r','LineStyle','--')
plot(lon_g(end,:),lat_g(end,:),'LineWidth',1.5,'Color','r','LineStyle','--')
plot(lon_g(:,1),lat_g(:,1),'LineWidth',1.5,'Color','r','LineStyle','--')
plot(lon_g(:,end),lat_g(:,end),'LineWidth',1.5,'Color','r','LineStyle','--')
plot(lon1089,lat1089,'LineWidth',1.5)
xlabel('Longitude')
ylabel('Latitude')
title([Casemean,' P1089'])
set(gca,'fontsize',14,'fontweight','b')

