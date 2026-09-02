function x1 = calcOracle(y,A,inputParam)
ind = find(inputParam.gg(:,1));
x1 = zeros(size(A,2),size(A,3));
A = A(:,ind,:);
numBS = size(A,3);
sigmax = 1/mean(inputParam.gammaPrior,'all');
% sigmaz = 1/inputParam.gammaOmega;
sigmaz = 0;

for tt = 1:numBS
    xtemp = (A(:,:,tt)'*A(:,:,tt) + sigmaz/sigmax*eye(length(ind)))*A(:,:,tt)'*y(:,tt);
    x1(ind,tt) = xtemp;
end
end