function H = defH(k, dk)
    % 定义从KE注入率误差到KE通量误差的变换矩阵
    % 输入:
    %   k: 波数数组
    %   dk: delta k数组
    % 输出:
    %   H: 变换矩阵
    
    nk = length(k);
    H = zeros(nk, nk+1);
    Hlog = zeros(nk, nk);
    
    for i = 1:nk
        Hlog(:, i) = double(k >= k(i)); % MATLAB中的heaviside等价
    end
    
    % 迭代所有元素
    for jj = 1:nk
        for ii = 1:nk
            H(jj, ii+1) = Hlog(jj, ii) * dk(jj);
        end
    end
    H(:, 1) = -ones(nk, 1);
return