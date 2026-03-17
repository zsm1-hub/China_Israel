function [R, Y, Y_std, ...
    optimal_fac, optimal_po, Th_Lag, Th_Lag_std, bootstrap_results] = ...
    Fk_fitting_SF3_Bayesian_RLS_Lcurve_uncertainty4(Case, nparticles, ...
    timerange, ini, mindist, maxdist, dot, n_bootstrap, bootstrap_size,kf1)
% input
% Case, nparticles, timerange, ini: for choose xxxx.mat with Lag SF3 and Lag CG
% dot: interval of kf (less points)
% mindist, maxdist: for choose range to fit
% n_bootstrap: bootstrap循环次数 (默认1000)
% bootstrap_size: 每次bootstrap采样样本数 (默认900)

% output
% R, Y, Y_std: Lag R, SF3, SF3_std       [dimension = r*1]
% kf, dkf: fitting point and diff
% Vt, eps, Fk: RLS-fitting SF3, energy injection rate, spectralfluxopen 
% filtscale, Th_lag, Th_lag_std: Coarse graining variables
% bootstrap_results: 包含所有bootstrap迭代结果的结构体

%% 1. 加载数据
fname = [Case, '_pars_P', num2str(nparticles), 'T', num2str(timerange(end)), ini, 'bootstrap.mat'];
load(fname);
Th_Lag=nanmean(Th_all,2);
Th_Lag_std=std(Th_all,0,2);
% 设置bootstrap参数
if nargin < 8
    n_bootstrap = 1000;
end
if nargin < 9
    bootstrap_size = 900;  % 从1000个样本中抽取900个
end

% 获取总样本数
n_total_samples = size(SF3, 2);  % SF3的维度应该是 r × n_samples

%% 2. 预计算几何参数（在循环外部）
% 先计算一次几何参数，用于确定数组大小
r = dist_axis;
ns = find(dist_axis >= mindist, 1);
ne = find(dist_axis <= maxdist, 1, 'last');

if isempty(ns) || isempty(ne)
    error('No valid points found in the specified distance range.');
end

R = r(ns:ne);
NR = length(R);

% 生成波数向量
kf_temp = logspace(log10(1/max(R)), log10(1/min(R)), dot).*2.*pi;
% kf_temp = logspace(log10(1/max(R)), log10(1/5000), dot).*2.*pi;

kf_temp=kf1;
% kf_temp = logspace(log10(1/max(R)), log10(1/min(R)), dot);

% kf_temp = logspace(log10(1/max(R)), log10(1/min(R))+log10(2.*pi), dot);
% kf_temp = logspace(log10(1/max(R)), log10(1/min(R)), dot);
% kf_temp = [logspace(log10(1/max(R)), log10(1/min(R)), dot),1./2e3.*2.*pi];

% kf_temp = logspace(log10(1/max(R)), log10(1/min(R)), dot);

% kf_temp = logspace(log10(1/max(R)), log10(1/min(R)), length(R)-1).*2.*pi;

% kf_temp = kf_temp(1:dot:end);
kf_temp = kf_temp';
dkf_temp = diff(kf_temp);
kf_final = 0.5 * (kf_temp(1:end-1) + kf_temp(2:end)); % 中点波数向量

% kf_final=kf1;
% 现在我们知道kf的大小，可以初始化数组了
n_kf = length(kf_final);
n_R = length(R);

% 初始化存储bootstrap结果的数组
bootstrap_results = struct();
bootstrap_results.eps_all = zeros(n_kf + 1, n_bootstrap);  % +1 用于epsilon_u
bootstrap_results.Fk_all = zeros(n_kf, n_bootstrap);
bootstrap_results.Vt_all = zeros(n_R, n_bootstrap);
% bootstrap_results.res_all = zeros(n_R, n_bootstrap);
bootstrap_results.optimal_fac_all = zeros(1, n_bootstrap);
bootstrap_results.optimal_po_all = zeros(1, n_bootstrap);
bootstrap_results.kf_final = zeros(n_kf,1);
% bootstrap_results.Fk_error_all = zeros(n_kf, n_bootstrap);
% bootstrap_results.eps_error_all = zeros(n_kf + 1, n_bootstrap);
% bootstrap_results.Cxx_all = cell(1, n_bootstrap);

%% 3. 主bootstrap循环
for kk = 1:n_bootstrap
    fprintf('Bootstrap iteration %d/%d\n', kk, n_bootstrap);
    
    try
        % 3.1 Bootstrap采样（有放回）
        bootstrap_indices = randi(n_total_samples, 1, bootstrap_size);
        SF3_bootstrap = SF3(:, bootstrap_indices);
        
        % 3.2 计算bootstrap样本的统计量
        SF3_Lag = nanmean(SF3_bootstrap, 2);
        SF3_Lag_std = std(SF3_bootstrap, 0, 2);
        
        % 3.3 提取当前距离范围内的数据
        Y = SF3_Lag(ns:ne);
        Y_std = std(SF3_bootstrap(ns:ne, :), 0, 2);
        
        % 3.4 使用预计算的几何参数
        R_use = R;
        kf_use = kf_final;
        dkf_use = dkf_temp;
        
        % 3.5 设置先验参数范围
        % fac0 = logspace(-17, -5, 80);  
        % po = logspace(-17, -5, 80); 
        % 
        % fac0 = logspace(-15, -5, 100);  
        % po = logspace(-15, -5, 100); 
        po =1e-7;
        fac0 = 1e-15;
        % 
 
        % 3.6 设置权重矩阵
        W = diag(Y_std.^2);
        
        % 3.7 构建设计矩阵A
        A = defA(R_use, kf_use, dkf_use);
        
        % 3.8 L曲线分析寻找最优先验参数
        [optimal_fac, optimal_po, results] = L_curve_analysis(Y, W, A, kf_use, fac0, po);
        
        % 存储最优参数
        bootstrap_results.optimal_fac_all(kk) = optimal_fac;
        bootstrap_results.optimal_po_all(kk) = optimal_po;
        
        % 3.9 设置先验协方差矩阵
        nk = length(kf_use);
        P_diag = [optimal_fac, ones(1, nk) * optimal_po];
        P = diag(P_diag);
        
        % 3.10 贝叶斯拟合
        [eps, Vt, res, Cxx] = RLS(Y, W, P, A);
        
        % 3.11 计算能谱通量
        Fk = calcFk(eps, kf_use, dkf_use);
        
        % 3.12 存储当前bootstrap迭代的结果
        bootstrap_results.eps_all(:, kk) = eps;
        bootstrap_results.Fk_all(:, kk) = Fk;
        bootstrap_results.Vt_all(:, kk) = Vt;
        bootstrap_results.kf_final=kf_use;
        % bootstrap_results.res_all(:, kk) = res;
        % bootstrap_results.Cxx_all{kk} = Cxx;
        
        % 3.13 计算不确定性
        % H = defH(kf_use, dkf_use);
        % Fxx = errorsFlux(Cxx, H);
        % bootstrap_results.Fk_error_all(:, kk) = sqrt(diag(Fxx));
        % bootstrap_results.eps_error_all(:, kk) = sqrt(diag(Cxx));
        
    catch ME
        fprintf('Bootstrap iteration %d failed: %s\n', kk, ME.message);
        % 存储NaN值表示失败
        bootstrap_results.eps_all(:, kk) = NaN;
        bootstrap_results.Fk_all(:, kk) = NaN;
        bootstrap_results.Vt_all(:, kk) = NaN;
        % bootstrap_results.res_all(:, kk) = NaN;
        bootstrap_results.optimal_fac_all(kk) = NaN;
        bootstrap_results.optimal_po_all(kk) = NaN;
        bootstrap_results.kf_final=NaN;
        % bootstrap_results.Fk_error_all(:, kk) = NaN;
        % bootstrap_results.eps_error_all(:, kk) = NaN;
        % bootstrap_results.Cxx_all{kk} = NaN;
    end
end
return