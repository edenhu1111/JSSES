function [dist_range, dist_vec] = calc_signal_dist_range(p_tx, p_rx, grid_matrix)
% 输入:
%   p_tx: 发射站位置 [x_s, y_s] (1 x 2)
%   p_rx: 接收站位置 [x_r, y_r] (1 x 2) 或 (N_rx x 2)
%   grid_matrix: 网格坐标矩阵 [x_g, y_g] (N_grid x 2)
% 输出:
%   dist_range: [min_dist, max_dist] 传播距离的全局范围
%   dist_vec: 每个网格点对应的双基距离向量 (N_grid x N_rx)

    N_grid = size(grid_matrix, 1);
    N_rx = size(p_rx, 1);
    dist_vec = zeros(N_grid, N_rx);

    % 计算发射站到所有网格的距离 (R1)
    d_tx_to_grid = sqrt(sum((grid_matrix - p_tx).^2, 2)); % (N_grid x 1)

    for i = 1:N_rx
        % 计算每个接收站到所有网格的距离 (R2)
        d_grid_to_rx = sqrt(sum((grid_matrix - p_rx(i,:)).^2, 2)); % (N_grid x 1)
        
        % 总路径距离 R = R1 + R2
        dist_vec(:, i) = d_tx_to_grid + d_grid_to_rx;
    end

    % 提取全局传播范围
    dist_range = [min(dist_vec(:))-10, max(dist_vec(:))+10];
end