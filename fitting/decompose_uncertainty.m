function [Cxx_epistemic, Cxx_aleatoric] = decompose_uncertainty(A, W, P, Cxx_total)
    % 分解参数协方差矩阵
    
    % 数据不确定性部分
    % Cxx_aleatoric = Cxx * A^T W^{-1} A * Cxx
    Winv = inv(W);
    M = A' * Winv * A;
    Cxx_aleatoric = Cxx_total * M * Cxx_total;
    
    % 认知不确定性部分  
    % Cxx_epistemic = Cxx * P^{-1} * Cxx
    Pinv = inv(P);
    Cxx_epistemic = Cxx_total * Pinv * Cxx_total;
    
    % 验证：Cxx_total ≈ Cxx_epistemic + Cxx_aleatoric
    % (在线性近似下成立)
end