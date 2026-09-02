function [pgOut] = mergeGrid(pgIn,x,xGrid,yGrid,d_min)
pg = pgIn.';

dist_sq = sum(pg.^2, 2) + sum(pg.^2, 2)' - 2*(pg * pg');
dist_sq(logical(eye(size(pg,1)))) = inf;

[row, col] = find(dist_sq < d_min^2);
ind = find(row < col);
row(ind) = [];
col(ind) = [];
for jj = 1:length(row)
    p_emerged(jj,:) = (pg(row(jj),:)*mean(abs(x(row(jj),:)).^2,2) + pg(col(jj),:)*mean(abs(x(col(jj),:)).^2,2))...
        /(mean(abs(x(row(jj),:)).^2,2) + mean(abs(x(col(jj),:)).^2,2));
    if norm(p_emerged(jj,:) - [xGrid(row(jj)),yGrid(row(jj))],2) < norm(p_emerged(jj,:) - [xGrid(col(jj)),yGrid(col(jj))],2)
        pg(row(jj),:) = p_emerged(jj,:);
        pg(col(jj),:) = [xGrid(col(jj)),yGrid(col(jj))];
    else
        pg(col(jj),:) = p_emerged(jj,:);
        pg(row(jj),:) = [xGrid(row(jj)),yGrid(row(jj))];
    end
end
pgOut = pg.';
end