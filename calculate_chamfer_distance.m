function cd_value = calculate_chamfer_distance(P1, P2)
% P1: 2 x N1 或 3 x N1 的第一组点云矩阵
% P2: 2 x N2 或 3 x N2 的第二组点云矩阵
% 返回值 cd_value: 对称的倒角距离（均方误差形式）

    % 确保输入维度正确 (坐标维度 x 点数)
    [dim1, N1] = size(P1);
    [dim2, N2] = size(P2);
    if dim1 ~= dim2
        error('两组点云的坐标维度必须一致（例如均为2D或3D）。');
    end

    % --- 计算第一项：P1 到 P2 的平均最近邻距离 ---
    % 采用分块或向量化计算以节省内存并保证速度
    % 对于每个 P1 中的点，寻找 P2 中最近的点
    dist_P1_to_P2 = compute_min_dist_sq(P1, P2);
    term1 = mean(dist_P1_to_P2);

    % --- 计算第二项：P2 到 P1 的平均最近邻距离 ---
    dist_P2_to_P1 = compute_min_dist_sq(P2, P1);
    term2 = mean(dist_P2_to_P1);

    % 对称 Chamfer Distance
    cd_value = term1 + term2;
end

function min_dist_sq = compute_min_dist_sq(A, B)
    % 计算 A 中每个点到 B 集合的最小距离平方
    % A: dim x Na, B: dim x Nb
    Na = size(A, 2);
    Nb = size(B, 2);
    
    % 如果数据量非常大，为防止内存溢出，采用循环分块处理
    % 如果数据量较小 (Na*Nb < 10^7)，可直接使用之前给您的矩阵化方法
    if Na * Nb < 5e6
        % 矩阵化计算所有点对的距离平方: |a-b|^2 = |a|^2 + |b|^2 - 2a'b
        D = sum(A.^2, 1)' + sum(B.^2, 1) - 2 * (A' * B);
        min_dist_sq = min(D, [], 2);
    else
        % 大规模数据分块处理
        min_dist_sq = zeros(Na, 1);
        blockSize = 1000; % 根据内存调整块大小
        for i = 1:blockSize:Na
            idx = i:min(i+blockSize-1, Na);
            A_block = A(:, idx);
            % 计算当前块 A 到全量 B 的距离
            D_block = sum(A_block.^2, 1)' + sum(B.^2, 1) - 2 * (A_block' * B);
            min_dist_sq(idx) = min(D_block, [], 2);
        end
    end
end