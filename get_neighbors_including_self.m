function unique_indices = get_neighbors_including_self(N, seed_indices)
    % N: 网格每边的点数 (N*N 网格)
    % seed_indices: 输入的一组一维序号 (1-indexed)
    
    % 1. 将一维序号转换为二维坐标 (r, c)
    % r: 行号 (1~N), c: 列号 (1~N)
    % 注意：MATLAB 的 ind2sub 默认也是列优先，可以直接使用
    [r, c] = ind2sub([N, N], seed_indices);
    
    % 2. 定义 8 邻域 + 自身 的相对偏移
    [dr, dc] = meshgrid(-1:1, -1:1);
    dr = dr(:); % 变为列向量
    dc = dc(:);
    
    % 3. 计算所有邻居的坐标
    % 利用广播机制：seed_r + dr'
    all_r = bsxfun(@plus, r(:), dr'); 
    all_c = bsxfun(@plus, c(:), dc');
    
    % 4. 展平并过滤越界坐标 (边界处理)
    all_r = all_r(:);
    all_c = all_c(:);
    
    valid_mask = (all_r >= 1 & all_r <= N & all_c >= 1 & all_c <= N);
    valid_r = all_r(valid_mask);
    valid_c = all_c(valid_mask);
    
    % 5. 将二维坐标转回一维序号
    neighbor_indices = sub2ind([N, N], valid_r, valid_c);
    
    % 6. 去重并排序
    unique_indices = unique(neighbor_indices);
end
