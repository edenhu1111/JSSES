function [RMSE, P_weighted, target_counts] = evaluate_offgrid_performance(P_meas, W_meas, P_true, d_min)
% P_meas: 2 x N_m 测量得到的坐标矩阵 (估计点)
% W_meas: 1 x N_m 测量点的权重 (通常是功率 |x_n|^2)
% P_true: 2 x N_t 真实的坐标矩阵 (Ground Truth)
% d_min:  判定阈值。距离真值小于 d_min 的点将被视为该真值的贡献分量
% -------------------------------------------------------------------------
% RMSE:   加权后的均方根误差
% P_weighted: 2 x N_t 每个真值点对应的加权测量中心 (若某真值无对应点则为 NaN)
% target_counts: 1 x N_t 每个真值点吸引到的测量点数量

    N_t = size(P_true, 2);
    N_m = size(P_meas, 2);
    
    % 初始化
    P_weighted = NaN(2, N_t);
    target_counts = zeros(1, N_t);
    sq_errors = []; % 存储有效匹配的平方误差

    % 1. 建立匹配关系 (Nearest Neighbor Association)
    % 对每一个测量点，寻找离它最近的真值点
    for i = 1:N_m
        % 计算该测量点到所有真值点的距离
        dists = sqrt(sum((P_true - P_meas(:, i)).^2, 1));
        [min_dist, nearest_idx] = min(dists);
        
        % 只有当距离小于阈值 d_min 时，才认为该测量值属于该真值
        if min_dist <= d_min
            % 记录该真值点“吸引”到的测量点信息
            % 我们先存储中间量，后续统一加权
            if target_counts(nearest_idx) == 0
                cluster_P{nearest_idx} = P_meas(:, i);
                cluster_W{nearest_idx} = W_meas(i);
            else
                cluster_P{nearest_idx} = [cluster_P{nearest_idx}, P_meas(:, i)];
                cluster_W{nearest_idx} = [cluster_W{nearest_idx}, W_meas(i)];
            end
            target_counts(nearest_idx) = target_counts(nearest_idx) + 1;
        end
    end

    % 2. 按照权重计算加权质心并统计误差
    for j = 1:N_t
        if target_counts(j) > 0
            % 执行加权中心计算: sum(P * W) / sum(W)
            weights = cluster_W{j};
            total_w = sum(weights);
            
            if total_w > 0
                P_weighted(:, j) = sum(cluster_P{j} .* weights, 2) / total_w;
            else
                P_weighted(:, j) = mean(cluster_P{j}, 2);
            end
            
            % 计算该加权中心与真值的平方误差
            sq_err = sum((P_weighted(:, j) - P_true(:, j)).^2);
            sq_errors = [sq_errors, sq_err];
        end
    end

    % 3. 计算 RMSE
    if isempty(sq_errors)
        RMSE = Inf;
        warning('没有测量点落入真值的 d_min 范围内，请检查阈值或算法收敛性。');
    else
        % RMSE = sqrt( mean( squared_errors ) )
        RMSE = sqrt(mean(sq_errors));
    end
end