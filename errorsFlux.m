function Fxx = errorsFlux(Cxx, H)
    % 从KE注入率误差转换为能量通量误差
    % 输入:
    %   Cxx: 包含KE注入率误差的方阵 (n_params × n_params)
    %   H: 变换矩阵 (n_k × n_params)
    % 输出:
    %   flux_std: 能量通量标准差向量 (n_k × 1)
    
    % 计算通量协方差矩阵
    Fxx = H * Cxx * H';
    
    % 提取对角线元素（方差）并开方得到标准差
    % flux_std = sqrt(diag(Fxx));
end