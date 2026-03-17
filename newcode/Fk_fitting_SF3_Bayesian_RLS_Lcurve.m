function [R,Y,Y_std,kf,dkf,Vt,eps,Fk,Fk_error,eps_error,...
    optimal_fac,optimal_po,...
    Th_Lag,Th_Lag_std]=Fk_fitting_SF3_Bayesian_RLS_Lcurve(Case,nparticles,...
    timerange,ini,mindist,maxdist,dot)
% input
% Case,nparticles,timerange,ini: for choose xxxx.mat with Lag SF3 and Lag CG
% dot: interval of kf (less points)
% mindist, maxdist: for choose range to fit

% output
% R, Y, Y_std: Lag R, SF3, SF3_std       [dimenison = r*1]
% kf, dkf: fitting point and diff
% Vt, eps, Fk: RLS-fitting SF3, energy injection rate, spectralflux
% filtscale, Th_lag, Th_lag_std: Coarse graining variables
%% 1. load variabels
fname=[Case,'_pars_P',num2str(nparticles),'T',num2str(timerange(end)),ini,...
    'bootstrap.mat'];
% load filtscale.mat
% filtscale=ncread(cgname,'filtscale');
% filtscale=filtscale(2:end);
load(fname);

Th_Lag=nanmean(Th_all,2);
Th_Lag_std=std(Th_all,0,2);
SF3_Lag=nanmean(SF3,2);
SF3_Lag_std=std(SF3, 0, 2);
SF3_Lag_all=(SF3);
r = dist_axis;
Nr = length(r);

% mindist=1e3;
% maxdist=500e3;
% 
ns = find(dist_axis >= mindist, 1);
ne = find(dist_axis <= maxdist, 1, 'last');

if isempty(ns) || isempty(ne)
    error('No valid points found in the specified distance range.');
end

R = r(ns:ne);
NR = length(R);
lf = R;


% if strcmp(kftype, 'log')
    kf = logspace(log10(1/max(R)), log10(1/min(R)), length(R)-1).*2.*pi;
% else
%     kf = kf1;
% end

kf=kf(1:dot:end);

kf = kf';
R=R';

Y=SF3_Lag(ns:ne);
Y_std=std(SF3(ns:ne,:), 0, 2);
% 
% fac0 = logspace(-15, -1, 80);  
% po = logspace(-15, -1, 80); 
fac0 = logspace(-15, -7, 80);  
po = logspace(-15, -7, 80); 
% fac0 = 1e-20;  
% po = 1e-20; 

dkf = diff(kf);
kf = 0.5*(kf(1:end-1) + kf(2:end)); % mid point wavenumber vec r*1

W = diag(Y_std.^2);
% W = diag((Y_std./Y_std).^2);

% 
A = defA(R, kf, dkf);

%% 2. L-curve to find best P choices
[optimal_fac, optimal_po, results] = L_curve_analysis(Y, W, A, kf, fac0, po);


nk = length(kf);
P_diag = [optimal_fac, ones(1, nk) * optimal_po];
% P_diag = [1e-5, ones(1, nk) * 3e-4];
P = diag(P_diag);

%% 3. fitting use Bayesian
[eps, Vt, res, Cxx] = RLS(Y, W, P, A);
Fk = calcFk(eps, kf, dkf);

% calc Fk errors and eps errors
H = defH(kf, dkf);
Fxx = errorsFlux(Cxx, H);
% FxxP = sqrt(diag(errorsFlux(P, H)));
% FxxD = sqrt(diag(errorsFlux(Cxx_alternative-P, H)));

Fk_error=sqrt(diag(Fxx));

eps_error=sqrt(diag(Cxx));

return