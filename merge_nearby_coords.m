function [P_merged, W_merged] = merge_nearby_coords(P, W, threshold)
% P: 2 x K 的矩阵，表示 K 个估计目标的 (x, y) 坐标
% W: 1 x K 的向量，表示每个目标对应的功率（或幅度的平方 |x_n|^2）
% threshold: 合并阈值，通常建议设置为 1.1 * grid_delta 到 1.5 * grid_delta
% ---------------------------------------------------------
% P_merged: 合并后的坐标矩阵
% W_merged: 合并后的总功率向量

    if isempty(P)
        P_merged = []; W_merged = []; return;
    end

    K = size(P, 2);
    visited = false(1, K);
    P_merged = [];
    W_merged = [];

    for i = 1:K
        if visited(i), continue; end
        
        % 1. 寻找与当前点距离在阈值内的所有点（包括自身）
        % 使用手动计算距离，避免依赖工具箱
        dists = sqrt(sum((P - P(:, i)).^2, 1));
        cluster_idx = find(dists <= threshold & ~visited);
        
        if ~isempty(cluster_idx)
            % 标记为已访问
            visited(cluster_idx) = true;
            
            % 2. 提取该簇的坐标和权重
            cluster_P = P(:, cluster_idx);
            cluster_W = W(cluster_idx);
            
            % 3. 执行功率加权质心计算 (Weighted Centroid)
            % 公式: P_new = sum(P_i * W_i) / sum(W_i)
            total_W = sum(cluster_W);
            
            if total_W > 0
                weighted_P = sum(cluster_P .* cluster_W, 2) / total_W;
            else
                % 若权重均为0，取算术平均
                weighted_P = mean(cluster_P, 2);
            end
            
            % 4. 存入结果
            P_merged = [P_merged, weighted_P];
            W_merged = [W_merged, total_W];
        end
    end
end