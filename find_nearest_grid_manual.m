function [nearest_grids, idx] = find_nearest_grid_manual(grid_points, target_points)
% grid_points: 2 x N
% target_points: 2 x K

    % 计算目标点每个分量的平方和 (1 x K)
    sum_target = sum(target_points.^2, 1); 
    % 计算网格点每个分量的平方和 (1 x N)
    sum_grid = sum(grid_points.^2, 1);     

    % 利用矩阵乘法计算中间项 -2*a'*b (K x N)
    % 这里计算的是每个目标点到每个网格点的距离平方
    % D(k, n) = |target_k|^2 + |grid_n|^2 - 2 * target_k' * grid_n
    D = sum_target' + sum_grid - 2 * (target_points' * grid_points);

    % 寻找每一行（即每个目标点）对应的最小距离索引
    [~, idx] = min(D, [], 2);
    idx = idx'; % 转为 1 x K

    % 提取结果
    nearest_grids = grid_points(:, idx);
end