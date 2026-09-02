function all_neighbors = get_neighbor_indices(ind, Nsize)
% 输入:
%   ind:   原始点的线性索引 (向量或标量)
%   Nsize: 网格的尺寸 (假设为 Nsize x Nsize)
% 输出:
%   all_neighbors: 包含原始点及其周围有效 8 邻域点的线性索引 (去重并排序)

    % 1. 将线性索引转换为二维坐标 [row, col]
    [r, c] = ind2sub([Nsize, Nsize], ind);
    
    % 2. 定义 8 邻域的相对偏移量 (包含中心点 [0,0] 则共 9 个点)
    [dr, dc] = meshgrid(-1:1, -1:1);
    dr = dr(:); dc = dc(:);
    
    neighbor_list = [];
    
    % 3. 遍历输入的每个索引，计算其邻居
    for i = 1:length(ind)
        % 计算潜在邻居的坐标
        rr = r(i) + dr;
        cc = c(i) + dc;
        
        % 4. 边界检查：确保坐标在 [1, Nsize] 范围内
        valid_mask = (rr >= 1 & rr <= Nsize) & (cc >= 1 & cc <= Nsize);
        valid_rr = rr(valid_mask);
        valid_cc = cc(valid_mask);
        
        % 5. 转换回线性索引
        neighbor_ind = sub2ind([Nsize, Nsize], valid_rr, valid_cc);
        neighbor_list = [neighbor_list; neighbor_ind]; %#ok<AGROW>
    end
    
    % 6. 去重并排序 (防止输入的点本身靠得很近导致索引重复)
    all_neighbors = unique(neighbor_list);
end