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
ini='_rough'
kftype='log'
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
% timerange=1:720;

% %500 m
% param_Eul_cg_sf3fk_500m
% SF3_E=SF3_mean';
% r_E=r;
% ymin_fk=-3e-8;ymax_fk=10e-8;
% ymin_sf=-3e-7;ymax_sf=1e-7;
% % ymin_ebsdk=-3e-8;ymax_ebsdk=3e-8;
% ymin_ebsdk=-2e-4;ymax_ebsdk=3e-4;
% ymin_right = -4e-8;    % 
% ymax_right = 6e-8;  % 
% xmin=0.5;
% xmax=1e3;

% % 2km param
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

% 2km param rot
% param_Eul_cg_sf3fk_rot
% SF3_E=nanmean(SF3(2:end,:),2);
% r_E=r(2:end)';
% ymin_fk=-3e-8;ymax_fk=3e-8;
% ymin_sf=-1e-7;ymax_sf=1e-7;
% % ymin_ebsdk=-3e-8;ymax_ebsdk=3e-8;
% ymin_ebsdk=-2e-4;ymax_ebsdk=3e-4;
% ymin_right = -4e-8;    % 
% ymax_right = 6e-8;  % 
% xmin=1;
% xmax=1e3;



if strcmpi(Case, 'wave')
    Casemean='Hf';
    if strcmpi(ini,'_rough')
        lambda=[1e-7,1e-8,1e-9];
        % str1=15;
        % end1=28;
        str1=1;
        end1=19;
        range1=2.5e3;
        range2=300e3;
    end
    if strcmpi(ini,'_roughsmall') || strcmpi(ini,'_roughsmall_rot') || strcmpi(ini,'_roughLASER') || strcmpi(ini,'_roughsmall_2mon') 
        % lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
        % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,1e-14];
        lambda=[1e-7,1e-8,1e-9];
        str1=1;
        end1=19;
        range1=2.5e3;
        range2=300e3;
    end
    if  strcmpi(ini,'_roughLASER') || strcmpi(ini,'_roughsmall_2mon') 
        % lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
        % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,1e-14];
        lambda=[1e-7,1e-8,1e-9];
        str1=1;
        end1=19;
        range1=2.5e3;
        range2=300e3;
    end
    if  strcmpi(ini,'_roughsmall_1mon') 
        % lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
        % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,1e-14];
        lambda=[1e-7,1e-8,1e-9];
        str1=1;
        end1=19;
        range1=2.5e3;
        range2=300e3;
        timerange=1:720;
    end
    if strcmpi(ini,'_roughsmall_500m')
        lambda=[1e-8];
        str1=1;
        end1=24;
        range1=2e3;
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
    if strcmpi(ini,'_roughsmall') || strcmpi(ini,'_roughsmall_rot') || strcmpi(ini,'_roughLASER') || strcmpi(ini,'_roughsmall_2mon') 
     %     lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
     % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,7e-14];
         lambda=[1e-7,1e-8,1e-9];
        str1=1;
        end1=18;
        range1=2.5e3;
        range2=300e3;
    end
    if strcmpi(ini,'_roughsmall_1mon') 
     %     lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
     % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,7e-14];
         lambda=[1e-7,1e-8,1e-9,1e-10,1e-11];
        str1=1;
        end1=18;
        range1=2.5e3;
        range2=300e3;
        timerange=1:720;
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

fname=[Case,'_pars_P',num2str(289),'T',num2str(timerange(end)),ini,'bootstrap.mat'];
load filtscale.mat

load(fname);
Th_Lag=nanmean(Th_all,2);
SF3_Lag=nanmean(SF3,2);
SF3_Lag_std=std(SF3, 0, 2);
SF3_Lag_all=(SF3);
r = dist_axis;
Nr = length(r);
mindist=2.5e3;
maxdist=500e3;
% 选择合理的距离范围
ns = find(dist_axis >= mindist, 1);
ne = find(dist_axis <= maxdist, 1, 'last');

if isempty(ns) || isempty(ne)
    error('No valid points found in the specified distance range.');
end

R = r(ns:ne);
NR = length(R);
lf = R;

% 创建波数向量
if strcmp(kftype, 'log')
    kf = logspace(log10(1/max(R)), log10(1/min(R)), length(R)-1).*2.*pi;
else
    kf = kf1;
end

% dk = diff(kf);
% kf = 0.5*(kf(1:end-1) + kf(2:end));
kf = kf';
% dk = dk';
R=R';
Y=SF3_Lag(ns:ne);

%%%%%%%%%%%%%%%%%%%%%%%%%fitting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Fk,k]=RLS_Bayesian(R,Y,kf,SF3_Lag_std(ns:ne))

semilogx(1./k.*2.*pi./1e3,Fk);hold on
semilogx(filtscale./1e3,Th_Lag);

semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{1});

semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end), ...
    'LineWidth',1.5, 'Color','k');

xlim([1,1e3])


%% L-curve
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
ini='_roughsmall'
kftype='log'
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
% timerange=1:720;

% %500 m
% param_Eul_cg_sf3fk_500m
% SF3_E=SF3_mean';
% r_E=r;
% ymin_fk=-3e-8;ymax_fk=10e-8;
% ymin_sf=-3e-7;ymax_sf=1e-7;
% % ymin_ebsdk=-3e-8;ymax_ebsdk=3e-8;
% ymin_ebsdk=-2e-4;ymax_ebsdk=3e-4;
% ymin_right = -4e-8;    % 
% ymax_right = 6e-8;  % 
% xmin=0.5;
% xmax=1e3;

% % 2km param
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

% 2km param rot
% param_Eul_cg_sf3fk_rot
% SF3_E=nanmean(SF3(2:end,:),2);
% r_E=r(2:end)';
% ymin_fk=-3e-8;ymax_fk=3e-8;
% ymin_sf=-1e-7;ymax_sf=1e-7;
% % ymin_ebsdk=-3e-8;ymax_ebsdk=3e-8;
% ymin_ebsdk=-2e-4;ymax_ebsdk=3e-4;
% ymin_right = -4e-8;    % 
% ymax_right = 6e-8;  % 
% xmin=1;
% xmax=1e3;



if strcmpi(Case, 'wave')
    Casemean='Hf';
    if strcmpi(ini,'_rough')
        lambda=[1e-7,1e-8,1e-9];
        % str1=15;
        % end1=28;
        str1=1;
        end1=19;
        range1=2.5e3;
        range2=300e3;
    end
    if strcmpi(ini,'_roughsmall') || strcmpi(ini,'_roughsmall_rot') || strcmpi(ini,'_roughLASER') || strcmpi(ini,'_roughsmall_2mon') 
        % lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
        % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,1e-14];
        lambda=[1e-7,1e-8,1e-9];
        str1=1;
        end1=19;
        range1=2.5e3;
        range2=300e3;
    end
    if  strcmpi(ini,'_roughLASER') || strcmpi(ini,'_roughsmall_2mon') 
        % lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
        % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,1e-14];
        lambda=[1e-7,1e-8,1e-9];
        str1=1;
        end1=19;
        range1=2.5e3;
        range2=300e3;
    end
    if  strcmpi(ini,'_roughsmall_1mon') 
        % lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
        % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,1e-14];
        lambda=[1e-7,1e-8,1e-9];
        str1=1;
        end1=19;
        range1=2.5e3;
        range2=300e3;
        timerange=1:720;
    end
    if strcmpi(ini,'_roughsmall_500m')
        lambda=[1e-8];
        str1=1;
        end1=24;
        range1=2e3;
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
    if strcmpi(ini,'_roughsmall') || strcmpi(ini,'_roughsmall_rot') || strcmpi(ini,'_roughLASER') || strcmpi(ini,'_roughsmall_2mon') 
     %     lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
     % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,7e-14];
         lambda=[1e-7,1e-8,1e-9];
        str1=1;
        end1=18;
        range1=2.5e3;
        range2=300e3;
    end
    if strcmpi(ini,'_roughsmall_1mon') 
     %     lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
     % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,7e-14];
         lambda=[1e-7,1e-8,1e-9,1e-10,1e-11];
        str1=1;
        end1=18;
        range1=2.5e3;
        range2=300e3;
        timerange=1:720;
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

[R,Y,Y_std,kf,dkf,Vt,eps,Fk,Fk_error,eps_error,...
    optimal_fac,optimal_po,filtscale,...
    Th_Lag,Th_Lag_std]=Fk_fitting_SF3_Bayesian_RLS_Lcurve(Case,1089,...
    timerange,ini,1e3,500e3,2);

% fname=[Case,'_pars_P',num2str(289),'T',num2str(timerange(end)),ini,'bootstrap.mat'];
% load filtscale.mat
% 
% load(fname);
% Th_Lag=nanmean(Th_all,2);
% Th_Lag_std=std(Th_all,0,2);
% SF3_Lag=nanmean(SF3,2);
% SF3_Lag_std=std(SF3, 0, 2);
% SF3_Lag_all=(SF3);
% r = dist_axis;
% Nr = length(r);
% mindist=1e3;
% maxdist=500e3;
% % 选择合理的距离范围
% ns = find(dist_axis >= mindist, 1);
% ne = find(dist_axis <= maxdist, 1, 'last');
% 
% if isempty(ns) || isempty(ne)
%     error('No valid points found in the specified distance range.');
% end
% 
% R = r(ns:ne);
% NR = length(R);
% lf = R;
% 
% % 创建波数向量
% if strcmp(kftype, 'log')
%     kf = logspace(log10(1/max(R)), log10(1/min(R)), length(R)-1).*2.*pi;
% else
%     kf = kf1;
% end
% kf=kf(1:2:end);
% 
% 
% kf = kf';
% R=R';
% 
% Y=SF3_Lag(ns:ne);
% Y_std=std(SF3(ns:ne,:), 0, 2);
% 
% % fac0 = logspace(-15, -1, 80);  
% % po = logspace(-15, -1, 80); 
% fac0 = logspace(-8, -3, 80);  
% po = logspace(-7, -4, 80); 
% 
% dkf = diff(kf);
% kf = 0.5*(kf(1:end-1) + kf(2:end)); % mid point wavenumber vec r*1
% 
% W = diag(Y_std.^2);
% 
% % 构建模型矩阵
% A = defA(R, kf, dkf);
% 
% [optimal_fac, optimal_po, results] = L_curve_analysis(Y, W, A, kf, fac0, po);
% 
% 
% nk = length(kf);
% P_diag = [optimal_fac, ones(1, nk) * optimal_po];
% P_diag = [1e-5, ones(1, nk) * 3e-4];
% 
% P = diag(P_diag);
% 
% [x0, y0, n0, Cxx] = RLS(Y, W, P, A);
% Fk = calcFk(x0, kf, dkf);
% 
% % 计算通量误差
% H = defH(kf, dkf);
% Fxx = errorsFlux(Cxx, H);
% % fk_err_laser_rls_mvb_logk = np.sqrt(np.diag(errorsFlux(cxx_laser_rls_mvb_logk, Herr_mvb_logk)))
% Fk_error=sqrt(diag(Fxx));
% 
subplot(2,2,1)
semilogx(1./kf.*2.*pi./1e3,Fk);hold on
semilogx(1./kf.*2.*pi./1e3,Fk+Fk_error);
semilogx(1./kf.*2.*pi./1e3,Fk-Fk_error);
semilogx(filtscale./1e3,Th_Lag);
% semilogx(filtscale./1e3,Th_Lag+Th_Lag_std);
% semilogx(filtscale./1e3,Th_Lag-Th_Lag_std);
semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{1});
semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end), ...
    'LineWidth',1.5, 'Color','k');
grid on
xlim([1,1e3])

subplot(2,2,2)
semilogx(R,Y./R);hold on
semilogx(R,Vt./R);
grid on

subplot(2,2,3)
% semilogx(1./kf.*(2*pi)./1e3, eps(2:end).*dkf)
semilogx(1./kf.*(2*pi)./1e3, eps(2:end))
grid on
xlim([1,1e3])

