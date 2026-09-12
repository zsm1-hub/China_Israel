function [x0, y0, res, Fk , k]=RLS_Bayesian(r,Y,k,Y_std,P)

dk = diff(k);
k = 0.5*(k(1:end-1) + k(2:end)); % mid point wavenumber vec r*1

W = diag(Y_std.^2);
% P = 1e-8 * eye(length(k)+1); % 先验协方差

% 构建模型矩阵
A = defA(r, k, dk);

% 执行正则化最小二乘
[x0, y0, res, Cxx] = RLS(Y, W, P, A);

% 计算能量通量
Fk = calcFk(x0, k, dk);

% 计算通量误差
H = defH(k, dk);
Fxx = errorsFlux(Cxx, H);
