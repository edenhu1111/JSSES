function [x1] = vampSVD(y,A,niter,inputParameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: vampSVD(BG prior)
% Description: VAMP(SVD form) for solving linear function 
% Author: Eden Hu (huyb@mail.sim.ac.cn)
%
% Input description
% y: the noisy data
% phi: sensing matrix
% gammaOmega: 
% niter: maximum iteration number
%
% Output description
% xest: estimated x
%               Nr->Na->Nv->Nsub
%               2D->1D: x outer, y inner (similar to vec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global x
gammaPrior = inputParameters.gammaPrior;
muPrior = inputParameters.muPrior;
Lambda = inputParameters.Lambda;
gammaOmega = inputParameters.gammaOmega;

M = size(A,1);
N = size(A,2);
% eta = @(x,beta) (x./abs(x)).*(abs(x)-beta).*(abs(x)-beta > 0); % denoising function
% dEta = @(x,beta) (abs(x)-beta > 0);
gamma1 = 0.01;

x1 = randn(N,1)/sqrt(gammaPrior);
% x1= x;
r1 = x1;
eps = 10^-16;
% dampFactorx = 1;
% dampFactoreta = 0.09;
[U,S,V] = svd(A); rA = rank(S);
U = U(:,1:rA);
S = S(1:rA,1:rA);
V = V(:,1:rA);
s = diag(S);
yTilde = diag(1./s)*U'*y;
for k = 1:niter
    xold = x1;
    gammaV = gamma1 + gammaPrior;
    m = (gamma1.*r1 + muPrior.*gammaPrior)./(gamma1 + gammaPrior);
    L = log(1./(1 + gamma1./gammaPrior)) + gamma1.*(r1.*conj(r1)) ...
    -((r1 - muPrior).*conj(r1 - muPrior)).*(gamma1.*gammaPrior)./(gamma1 + gammaPrior);
    piL = Lambda./(Lambda + (1-Lambda).*exp(-L));
    x1 = piL.*m;
    eta1 = 1./mean(piL.*(m.*conj(m) + 1/gammaV)-piL.^2.*(x1.*conj(x1)));
    gamma2 = max(mean(eta1 - gamma1),eps);
    r2 = (eta1.*x1 - gamma1.*r1)./gamma2;

    
    d = gammaOmega*diag(1./(gammaOmega*s.^2 + gamma2*ones(rA,1)))*s.^2;
    gamma1 = gamma2*sum(d)./(N - sum(d));
    r1 = r2 + N*V*diag(d ./ sum(d) )*(yTilde - V'*r2);
    if norm(x1-xold,'fro') < 1e-16
        break;
    end
end
end