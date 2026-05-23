function clean_pairs_time=clear_nan_in_pairs_time(pairs_time)
clean_pairs_time = struct('dist', {}, 'dul', {}, 'dut', {},'theta',{});

% 遍历每个时间点
for i = 1:length(pairs_time)
    % 获取当前时间点的数据
    current_dist = pairs_time(i).dist;
    current_dul = pairs_time(i).dul;
    current_dut = pairs_time(i).dut;
    current_theta = pairs_time(i).theta;
    
    % 找出所有字段都非NaN的位置
    valid_idx = ~isnan(current_dist) & ...
                ~isnan(current_dul) & ...
                ~isnan(current_dut) & ...
                ~isnan(current_theta);
    
    % 只保留完全有效的点
    clean_dist = current_dist(valid_idx);
    clean_dul = current_dul(valid_idx);
    clean_dut = current_dut(valid_idx);
    clean_theta = current_theta(valid_idx);
    
    % 将处理后的数据存入新结构体
    clean_pairs_time(i).dist = clean_dist;
    clean_pairs_time(i).dul = clean_dul;
    clean_pairs_time(i).dut = clean_dut;
    clean_pairs_time(i).theta = clean_theta;
    
    % 显示进度（可选）
    if mod(i, 100) == 0
        fprintf('已处理 %d/%d 个时间点\n', i, length(pairs_time));
    end
end
return
