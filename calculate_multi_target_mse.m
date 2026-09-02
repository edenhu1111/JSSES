function [mse, matched_est, mappings] = calculate_multi_target_mse(est_pos,true_pos)
    % 输入:
    %   true_pos: 真值坐标矩阵 [K x 2]
    %   est_pos:  算法估计的坐标池 [M x 2]
    % 输出:
    %   mse:         最优匹配后的均方误差
    %   matched_est: 与真值一一对应的估计点 [K x 2]
    %   mappings:    匹配索引 [K x 2]，第一列为真值索引，第二列为估计值索引

    K = size(true_pos, 1);
    M = size(est_pos, 1);
    
    % 1. 计算所有点对之间的欧氏距离平方 (手动实现，无需 pdist2)
    % 利用广播机制: (K x 1 x 2) - (1 x M x 2)
    diff = reshape(true_pos, K, 1, []) - reshape(est_pos, 1, M, []);
    dist_sq_matrix = sum(diff.^2, 3); 

    % 2. 执行最优二分匹配 (匈牙利算法)
    % 如果没有工具箱，这里使用 matchpairs。若无此函数，可替换为下方的简易版本。
    try
        % matchpairs 寻找总距离平方和最小的匹配方式
        idx_pairs = matchpairs(dist_sq_matrix, 1e6); 
    catch
        % 如果环境没有 matchpairs，采用贪心匹配（虽非全局最优，但对高 SNR 场景足够）
        idx_pairs = greedy_match_logic(dist_sq_matrix);
    end

    % 3. 提取匹配点并计算 MSE
    matched_true_idx = idx_pairs(:, 1);
    matched_est_idx = idx_pairs(:, 2);
    
    matched_est = est_pos(matched_est_idx, :);
    actual_true = true_pos(matched_true_idx, :);
    
    % 计算 MSE = (1/K) * sum( ||true - est||^2 )
    errors_sq = sum((actual_true - matched_est).^2, 2);
    mse = sum(errors_sq) / size(idx_pairs, 1);
    mappings = idx_pairs;

    % 检查是否有真值未被匹配
    if size(idx_pairs, 1) < K
        warning('部分真值未找到匹配的估计点！');
    end
end

function idx_pairs = greedy_match_logic(cost_mat)
    % 简易贪心匹配：逐个找出全局最小距离并剔除行列
    [K, M] = size(cost_mat);
    idx_pairs = [];
    temp_mat = cost_mat;
    for i = 1:min(K, M)
        [min_val, min_idx] = min(temp_mat(:));
        if isinf(min_val), break; end
        [r, c] = ind2sub(size(temp_mat), min_idx);
        idx_pairs = [idx_pairs; r, c];
        temp_mat(r, :) = inf; % 剔除已匹配的真值
        temp_mat(:, c) = inf; % 剔除已匹配的估计点
    end
end