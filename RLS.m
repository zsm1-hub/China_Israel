function [x0, y0, n0, Cxx] = RLS(Y, W, P, X)
    % 正则化最小二乘
    % 输入:
    %   Y = 数据
    %   W = 权重矩阵（2D或1D数组）
    %   P = 参数不确定性（2D或1D数组）
    %   X = 要拟合的模型
    %
    % 输出:
    %   x0 = 系数
    %   y0 = 拟合值
    %   n0 = 残差
    %   Cxx = 不确定性

    % 检查W矩阵维度
    if ndims(W) == 2
        if size(W, 1) ~= size(W, 2)
            error('W不是方阵');
        elseif size(diag(W), 1) ~= size(X, 1)
            error('diag(W)长度与X行数不匹配');
        end
    elseif isvector(W)
        W = diag(W);
        if size(W, 1) ~= size(X, 1)
            error('diag(W)长度与X行数不匹配');
        end
    end
    
    % 检查P矩阵维度
    if isvector(P)
        if length(P) == 1
            P = P * eye(size(X, 2));
        else
            if length(P) == size(X, 2)
                P = diag(P);
            else
                error('diag(P)长度与X列数不匹配');
            end
        end
    elseif ndims(P) == 2
        if size(P, 1) ~= size(P, 2)
            error('P不是方阵');
        elseif length(diag(P)) ~= size(X, 2)
            error('diag(P)长度与X列数不匹配');
        end
    end
    
    % 检查Y和X维度匹配
    if length(Y) ~= size(X, 1)
        error('Y长度与X行数不匹配');
    end
    
    % 计算逆矩阵
    if mean(P(:)) == 0
        Pinv = P;
    else
        Pinv = inv(P); % P^{-1}
    end
    
    Winv = inv(W); % W^{-1}
    
    % 计算正则化解
    M = X' * Winv * X; % M = (A' x W^{-1} x A)
    Minv = inv(M + Pinv); % Minv = (M + P^{-1})^{-1}
    N = Minv * X' * Winv; % N = Minv x A' x W^{-1}
    
    % 解
    x0 = N * Y; % x0 = N x Y
    
    % 残差
    y0 = X * x0;
    n0 = Y - y0; % n0 = Y - A x x0
    
    % 不确定性
    Cxx = Minv;
return