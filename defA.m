function A = defA(r, k, dk)
    % 定义模型矩阵
    % 输入:
    %   r: 距离分箱
    %   k: 波数分箱
    %   dk: delta k
    % 输出:
    %   A: 模型矩阵
    
    nr = length(r);
    nk = length(k);
    A = zeros(nr, nk+1);
    
    % 迭代所有元素
    for jj = 1:nk
        for ii = 1:nr
            A(ii, jj+1) = -4 * besselj(1, k(jj)*r(ii)) / k(jj) * dk(jj);
        end
    end
    
    A(:, 1) = 2 * r;
return