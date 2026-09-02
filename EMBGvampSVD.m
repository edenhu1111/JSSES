function [x1,state] = EMBGvampSVD(y,A,inputParam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: EMBGvampSVD(BG prior)
% Description: VAMP(SVD form) for solving linear function 
% Author: Eden Hu (huyb@mail.sim.ac.cn)
%
% Input description
% y: the noisy data
% phi: sensing matrix 
% inputParam:
% -gammaOmega:inverse of noise power
% -niter: maximum iteration number 
% -EMiter:maximum EM iteration number
% -muPrior: The initial mean value of BG prior 
% -Lambda: initial bernoulli parameter 
% -gammaPrior: initial variation of BG prior
%
% Output description
% xest: estimated x 
% state: 
% -NMSE: the NMSE of the estimated x 
% -sparsity: the average sparsity of x
%               Nr->Na->Nv->Nsub
%               2D->1D: x outer, y inner (similar to vec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global gg
gammaPrior = mean(inputParam.gammaPrior,'all');
muPrior = mean(inputParam.muPrior,'all');
Lambda = inputParam.Lambda;
gammaOmega = inputParam.gammaOmega;
niter = inputParam.niter;
EMiter = inputParam.EMiter;
convBreaker = inputParam.convBreaker;
dampFacGam = inputParam.dampFacGam;
dampFac = inputParam.dampFac;

gg = inputParam.gg;

M = size(A,1);
N = size(A,2);
NMSE = zeros(EMiter,1);
NMSEz = zeros(niter,EMiter);

sparsity = zeros(niter,EMiter);

gamma2 = mean(gammaPrior);


x1 = Lambda*muPrior + (randn(N,1)+ 1j*randn(N,1))/sqrt(gammaPrior/2);
% r1 = zeros(N,1);
r2 = zeros(N,1);
eps = 10^-16;
[U,S,V] = svd(A); 
% rA = rank(S);
ss=diag(S); rA = sum((ss/max(ss))>=1e-5);
U = U(:,1:rA);
S = S(1:rA,1:rA);
V = V(:,1:rA);
s = diag(S);
yTilde = diag(1./s)*U'*y;
% AA = A.*conj(A);
gammaMax = 1e+15;
piL = Lambda*ones(N,1);
for t = 1:EMiter
%% E step    
    xoldOut = x1;
    for k = 1:niter
        %% Equivalent LMMSE Estimator(EXT message for MMSE estimator)
%         cov2 = inv(gammaOmega*(A'*A) + gamma2*eye(N));
%         x2 = cov2*(gammaOmega*A'*y + gamma2*r2);
%         alpha2 = gamma2/N*trace(cov2);
%         eta2 = min(max(gamma2/alpha2,1/gammaMax),gammaMax);
%         gamma1 = eta2 - gamma2;
%         r1 = (eta2*x2 - gamma2.*r2)./gamma1;

        d = ( gammaOmega*s.^2./(gammaOmega*s.^2 + gamma2*ones(rA,1)) );

        x2 = r2 + V*diag(d)*(yTilde - V'*r2);
        muz = sum( s.^2 ./ (gammaOmega*s.^2 + gamma2))/M;

        gamma1 = gamma2*sum(d)./(N - sum(d));
        r1 = r2 + (N/rA)*V*diag(d/mean(d))*(yTilde - V'*r2);

        gamma1 = max(min(gammaMax,gamma1),1/gammaMax);
        % 
        
        %% Damping
        if k > 1 && t > 1 
            r1 = dampFac*r1 + (1-dampFac)*r1Old;
            gamma1 = dampFacGam*(gamma1) + (1-dampFacGam)*gamma1old;
        elseif t > 1 && k == 1
            gamma1 = dampFacGam*(gamma1) + (1-dampFacGam)*gamma1old;
        end
        r1Old = r1; gamma1old = gamma1;

        %% MMSE estimator
        
        xold = x1;
        gammaV = gamma1 + gammaPrior;
        m = (gamma1.*r1   + muPrior.*gammaPrior)./(gamma1 + gammaPrior);
        % m(isinf(gammaPrior)) = 0;m(isnan(gammaPrior)) = 0;
        
        sigmaSum = 1./gammaPrior + 1./gamma1;
        NNum  = Lambda./sigmaSum/pi;
        NDeno = (1-Lambda).*exp(abs(r1 - muPrior).^2./sigmaSum-abs(r1).^2.*gamma1).*gamma1/pi;
        piL = NNum./(NNum + NDeno);
        x1 = piL.*m;
        eta1 = 1./mean( piL.*(abs(m).^2 + 1./gammaV) - abs(x1).^2 );

        %% EXT message from MMSE to LMMSE
        % gamma2 = max(mean(eta1 - gamma1),eps);
        gamma2 = eta1 - gamma1;
        gamma2 = max(min(gammaMax,gamma2),1/gammaMax);
        r2 = x1 + gamma1.*(x1 - r1)./(eta1 - gamma1);

        %% Damping
        if t > 1 && k > 1
            r2 = dampFac*r2 + (1-dampFac)*r2Old;
            gamma2 = dampFacGam*gamma2 + (1-dampFacGam)*gamma2Old;
        end
        gamma2Old = gamma2; r2Old = r2;

        %% record convergence information


        sparsity(k,t) = mean(piL); 
        NMSEz(k,t) = norm(y - A*x1,2)^2/norm(A*gg,2)^2;
        state.eta1(k,t) = eta1;
        state.gamma1(k,t) = gamma1;
        state.gamma2(k,t) = gamma2;
        if norm(x1-xold,'fro')^2/norm(xold,'fro')^2 < 1e-8 && convBreaker
            break;
        end
    end
 %% M step
    gammaOmega = 1./(mean((y - A*x2).*conj(y - A*x2)) + muz);
%     Lambda = max(mean(piL),0); 
%     gammaPrior = 1./( ( (muPrior - m).*conj(muPrior - m) + 1./gammaV ) );
%     muPrior = m;
    
    %%%% EM Maximization
    Lambda = max(mean(piL),0);
    gammaPrior = Lambda./mean(piL.*( (muPrior - m).*conj(muPrior - m) + 1./gammaV));
%     muPrior = mean(piL.*m)/Lambda;
    if inputParam.Normalization == 1
        xi = (x1'*gg)/(gg'*gg);
    else
        xi = 1;
    end
    NMSE (t) = norm(x1 - xi*gg,2)^2/norm(gg,2)^2;
%     muPrior = m;
    if norm(x1 - xoldOut,'fro')/norm(xoldOut,'fro') < 1e-8 && convBreaker
        break;
    end
    %% record EM's M-step result
%     state.mu(t) = norm(muPrior,'fro');
    state.gammaP(t) = mean(Lambda./gammaPrior);
    state.lambda(t) = mean(Lambda);
    state.gammaOmega(t) = gammaOmega;
    state.muP(:,t) = muPrior;
end
state.NMSE = NMSE;
state.NMSEz = NMSEz;
state.sparsity = sparsity;
end