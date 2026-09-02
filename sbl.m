function [x, state] = sbl(y,Phi, inputParam)
%SBL Function to apply the sbl (EM) algorithm
%   y: measurements, Phi: measurement matrix, sigma: noise std. The
%   algorithm is executed till ||y - Phi x|| < eps if Nmax = 0 else Nmax iterations. Returns sparse x and
%   time elapsed
tic
global gg
%% Initialization
M = size(Phi,1);
N = size(Phi,2);

eps = 1e-6;

rho0 = inputParam.rho0;
Niter = inputParam.Niter;
a = inputParam.gammaA;
b = inputParam.gammaB;

c = inputParam.gammaC;
d = inputParam.gammaD;

rho = ones(N,1);
x = zeros(N,1);

[~,singV,~] = svd(Phi);
singV = diag(singV);
%% EM Iteration
for niter = 1:Niter
    xOld = x;
    % update X
    S = inv((Phi'*Phi)*rho0 + diag(rho));
    x = S*(Phi'*y)*rho0;

%     muz = sum(singV.^2 ./ (rho0.*singV.^2 + rho))/M;
    muz = trace(Phi*S*Phi')/M;
    ss = abs(diag(S));
    %update rho and rho0
    rho =  (a+1)./(abs(x).^2 + ss + b);
    rho0 = (c+1)./(norm(y - Phi*x,2)^2/M + muz +d);
    
    state.NMSE(niter) = norm(x - gg,2)^2/norm(gg)^2;

    % Break when converging
    if norm(x - xOld,2) < eps
        break;
    end
end
state.rho = rho;
state.rho0 = rho0;
telap = toc;
end