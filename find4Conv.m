function [x] = find4Conv(X)
N = size(X,2);
ind = zeros(N,1);
for nn = 1:N
    ind(nn) = find(X(:,nn),1,'last');
    x(nn) = X(ind(nn),nn);
end

end