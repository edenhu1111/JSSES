function [x1] = vamp(y,A,niter,gammaPrior,muPrior,Lambda,gammaOmega)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: vamp(laplacian prior)
% Description: VAMP for solving linear function 
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
M = size(A,1);
N = size(A,2);
% eta = @(x,beta) (x./abs(x)).*(abs(x)-beta).*(abs(x)-beta > 0); % denoising function
% dEta = @(x,beta) (abs(x)-beta > 0);
gamma1 = 0.01;

% gammaPrior = 0.01;
% muPrior = 0;
% Lambda = 0.05;
% r1 = (randn(N,1)+1j*randn(N,1))/sqrt(2);
% r1 = sqrt(1/gamma1)*randn(N,1);
x1 = randn(N,1)/sqrt(gammaPrior);
% x1= x;
r1 = x1;
eps = 10^-16;
% dampFactorx = 1;
% dampFactoreta = 0.09;

for k = 1:niter
    % x1 = dampFactorx*eta(r1,gamma1) + (1-dampFactorx)*x1;
    % alpha1 = dampFactor*mean(dEta(r1,gamma1)) + (1-dampFactor)*alpha1;
    % alpha1 = mean(dEta(r1,gamma1));
    xold = x1;
    % x1 = (gamma1*r1)/(gamma1 + gammaPrior)*(1-Lambda);
    % eta1 =  (gamma1 + gammaPrior)/(1-Lambda) ;
    gammaV = gamma1 + gammaPrior;
    m = (gamma1.*r1 + muPrior.*gammaPrior)./(gamma1 + gammaPrior);
    % L = 1/2*log(1./(1 + gamma1./gammaPrior)) + gamma1.*(r1.*conj(r1))./2 ...
    %     -1/2.*(r1.*conj(r1)).*(gamma1*gammaPrior)./(gamma1+gammaPrior);
    L = log(1./(1 + gamma1./gammaPrior)) + gamma1.*(r1.*conj(r1)) ...
    -((r1 - muPrior).*conj(r1 - muPrior)).*(gamma1.*gammaPrior)./(gamma1 + gammaPrior);
    piL = Lambda./(Lambda + (1-Lambda).*exp(-L));
    x1 = piL.*m;
    eta1 = 1./mean(piL.*(m.*conj(m) + 1/gammaV)-piL.^2.*(x1.*conj(x1)));
    % eta1 = max(gamma1/alpha1,eps);
    % eta1 = dampFactoreta*gamma1/alpha1 + (1-dampFactoreta)*eta1;
    gamma2 = max(mean(eta1 - gamma1),eps);
    r2 = (eta1.*x1 - gamma1.*r1)./gamma2;

    
    % cov2 = inv(gammaOmega*(A'*A) + gamma2*eye(N));
    cov2 = 1/gamma2*eye(N) - ...
        gammaOmega/gamma2^2*A'*inv(eye(M)+gammaOmega/gamma2*(A*A'))*A;       % Woodbury formula
    x2 = cov2*(gammaOmega*A'*y + gamma2*r2);
    alpha2 = gamma2/N*trace(cov2);
    eta2 = max(gamma2/alpha2,eps);
    gamma1 = max(eta2 - gamma2,eps);
    r1 = (eta2*x2 - gamma2.*r2)./gamma1;
    if norm(x1-xold,'fro') < 1e-16
        break;
    end
end
end