function saveResults_pairs(fname, ii, nvaild, ...
                    pairs_time, total_steps)
    % 创建文件名
    save_fname = [fname(1:end-3), 'chunk', num2str(ii, '%02d'), 'pairs.mat'];
    p1=pairs_time;
    clear pairs_time;
    pairs_time=clear_nan_in_pairs_time(p1);
    clear p1

    % 保存数据
    save(save_fname, 'nvaild', ...
         "pairs_time", 'total_steps', '-v7.3');
end