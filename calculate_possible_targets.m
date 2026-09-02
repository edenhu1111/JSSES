function targets = calculate_possible_targets(p1, phi1, theta1, p2, phi2, theta2)
    % 输入:
    % p1, p2: 车站位置 [x; y]
    % phi1, phi2: 天线偏置角 (弧度)
    % theta1: 站1观测到的AoA数组 (1 x K1)
    % theta2: 站2观测到的AoA数组 (1 x K2)
    % 输出:
    % targets: 2 x (K1*K2) 的矩阵，包含所有可能的组合坐标
    
    K1 = length(theta1);
    K2 = length(theta2);
    
    % 预分配空间：总共有 K1 * K2 种组合
    targets = NaN(2, K1 * K2);
    
    % 转换到全局坐标系的绝对方位角
    alpha1 = theta1 + phi1;
    alpha2 = theta2 + phi2;
    
    % 预计算方向向量，提高循环效率
    U1 = [cos(alpha1); sin(alpha1)]; % 2 x K1
    U2 = [cos(alpha2); sin(alpha2)]; % 2 x K2
    
    b = p2 - p1;
    
    for i = 1:K1
        u1 = U1(:, i);
        for j = 1:K2
            u2 = U2(:, j);
            
            % 求解方程组: d1*u1 - d2*u2 = p2 - p1
            % [u1, -u2] * [d1; d2]' = b
            A = [u1, -u2];
            
            % 检查奇异性（平行线）
            if abs(det(A)) > 1e-6
                dist = A \ b;
                % 物理可行性判断：d1, d2 必须为正（目标在射线前方）
                if dist(1) > 0 && dist(2) > 0
                    targets(:, (i-1)*K2 + j) = p1 + dist(1) * u1;
                end
            end
        end
    end
end