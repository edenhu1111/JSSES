function [centroids, cluster_energies, count] = cluster_weighted_centroid(pts_xy, weights, epsilon)
% 输入:
%   pts_xy:  坐标矩阵 (N x 2)
%   weights: 权重向量 (N x 1)
%   epsilon: 距离阈值 (建议取网格步长的 1.1~1.5 倍)
%
% 输出:
%   centroids:        加权质心 (K x 2)
%   cluster_energies: 每个簇的总能量 (K x 1)
%   count:            目标个数
    % 1. 预处理：只保留有能量的点
    th = max(weights) * 0.01; % 剔除能量极低的点
    keep = weights > th;
    pts = pts_xy(keep, :);
    w = weights(keep);
    
    if isempty(pts)
        centroids = []; cluster_energies = []; count = 0; return;
    end

    N = size(pts, 1);
    visited = false(N, 1);
    cluster_id = zeros(N, 1);
    curr_id = 0;

    % 2. 简易聚类逻辑 (类似 BFS 搜索)
    for i = 1:N
        if ~visited(i)
            curr_id = curr_id + 1;
            % 寻找与当前点距离小于 epsilon 的所有点
            queue = i;
            visited(i) = true;
            
            head = 1;
            while head <= length(queue)
                curr_node = queue(head);
                cluster_id(curr_node) = curr_id;
                
                % 计算当前点到其他所有未访问点的距离
                dist_sq = sum((pts - pts(curr_node,:)).^2, 2);
                neighbors = find(dist_sq < epsilon^2 & ~visited);
                
                % 入队并标记
                visited(neighbors) = true;
                queue = [queue; neighbors]; %#ok<AGROW>
                head = head + 1;
            end
        end
    end

    % 3. 计算加权质心
    count = curr_id;
    centroids = zeros(count, 2);
    cluster_energies = zeros(count, 1);

    for k = 1:count
        idx_k = (cluster_id == k);
        w_k = w(idx_k);
        pts_k = pts(idx_k, :);
        
        sum_w = sum(w_k);
        % 加权质心公式
        centroids(k, :) = sum(w_k .* pts_k, 1) / (sum_w);
        cluster_energies(k) = sum_w;
    end
    
    % 按能量从大到小排序
    [cluster_energies, s_idx] = sort(cluster_energies, 'descend');
    centroids = centroids(s_idx, :);
end