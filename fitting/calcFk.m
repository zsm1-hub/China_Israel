function Fk = calcFk(eps, k, dk)
    % 计算动能通量
    % 输入:
    %   eps: KE注入率向量 (n_params × 1)，第一个元素是eps_u，后续是eps_j
    %   k: 波数向量 (n_k × 1)
    %   dk: 波数间隔向量 (n_k × 1)
    % 输出:
    %   Fk: 动能通量向量 (n_k × 1)
    
    nk = length(k);
    Fk = zeros(nk, 1);  % 确保是列向量
    
    % 第一个波数的通量
    Fk(1) = -eps(1);
    
    % 累加后续波数的贡献
    for jj = 1:nk-1
        Fk(jj+1) = Fk(jj) + eps(jj+1) * dk(jj);
    end
end