function [Fk_total_std, Fk_epistemic_std, Fk_aleatoric_std] = decompose_Fk_uncertainty(Cxx, H, A, W, P)
    % Fk谱的不确定性分解
    
    % 1. 首先分解参数不确定性
    [Cxx_epistemic, Cxx_aleatoric] = decompose_uncertainty(A, W, P, Cxx);
    
    % 2. 传播到Fk谱
    % 总不确定性
    C_Fk_total = H * Cxx * H';
    Fk_total_std = sqrt(diag(C_Fk_total));
    
    % 认知不确定性
    C_Fk_epistemic = H * Cxx_epistemic * H';
    Fk_epistemic_std = sqrt(diag(C_Fk_epistemic));
    
    % 数据不确定性  
    C_Fk_aleatoric = H * Cxx_aleatoric * H';
    Fk_aleatoric_std = sqrt(diag(C_Fk_aleatoric));
    
    % 验证：总方差 ≈ 认知方差 + 数据方差
    % (在协方差矩阵可加的条件下)
end