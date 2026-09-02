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
global gg fc c0
gammaPrior = inputParam.gammaPrior;
muPrior = inputParam.muPrior;
Lambda = inputParam.Lambda;
gammaOmega = inputParam.gammaOmega;
niter = inputParam.niter;
EMiter = inputParam.EMiter;
EMOutIter = inputParam.EMOutiter;
convBreaker = inputParam.convBreaker;
dampFacGam = inputParam.dampFacGam;
dampFac = inputParam.dampFac;
subcarrInd = inputParam.subcarrInd;


xUAV  = inputParam.xUAV;

sigmaPos = inputParam.sigmaPos;
Tdelay  = inputParam.Tdelay;
xR = inputParam.xR;
stepSize = inputParam.stepSize;
pilot = ones(length(subcarrInd),1);

armijoBeta = inputParam.armijoBeta;
armijoSigma = inputParam.armijoSigma;

Nr = size(xR,2);
M = size(A,1);
N = size(A,2);

Hdx = inputParam.sMDeriv.Phidx(:,1:N);
Hdy = inputParam.sMDeriv.Phidy(:,1:N);
Hdz = inputParam.sMDeriv.Phidz(:,1:N);


NMSE = zeros(niter,EMiter);
NMSEz = zeros(niter,EMiter);
% state.mu = zeros(EMiter,1);
% state.gammaP = zeros(EMiter,1);
% state.lambda = zeros(EMiter,1);
sparsity = zeros(niter,EMiter);
% state.eta1 = zeros(niter,EMiter);
% state.gamma1 = zeros(niter,EMiter);
% state.gamma2 = zeros(niter,EMiter);
% state.gammaOmega = zeros(EMiter,1);
% state.xRec = zeros(N,3*niter);
% % state.mRec = zeros(N,3*niter);
% state.xrRec = zeros(N,3*niter);
% state.muP = zeros(N,EMiter);
% gamma1 = gammaPrior;
gamma2 = gammaPrior;

gammaMax = 1e+15;
piL = Lambda*ones(N,1);

% r1 = zeros(N,1);
for tt = 1:EMOutIter
    x1 = Lambda*muPrior + (randn(N,1)+ 1j*randn(N,1))/sqrt(gammaPrior/2);
    r2 = zeros(size(x1));
%     eps = 10^-16;
    [U,S,V] = svd(A); 
    % rA = rank(S);
    ss=diag(S); rA = sum(ss>=0.001);
    U = U(:,1:rA);
    S = S(1:rA,1:rA);
    V = V(:,1:rA);
    s = diag(S);
    yTilde = diag(1./s)*U'*y;
    % AA = A.*conj(A);
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
            end
            r1Old = r1; gamma1old = gamma1;
    
            %% MMSE estimator
            
            xold = x1;
            gammaV = gamma1 + gammaPrior;
            m = (gamma1.*r1   + muPrior.*gammaPrior)./(gamma1 + gammaPrior);
            % m(isinf(gammaPrior)) = 0;m(isnan(gammaPrior)) = 0;
            
            sigmaSum = 1./gammaPrior + 1./gamma1;
            NNum  = Lambda.*exp(-abs(r1 - muPrior).^2./sigmaSum)./sigmaSum/pi;
            NDeno = (1-Lambda).*exp(-abs(r1).^2.*gamma1).*gamma1/pi;
            piL = NNum./(NNum + NDeno);
    %         Lmax = 1e5;
    %         L = max(min(L,Lmax),-Lmax);
    %         piL = 1./(1 + (1-Lambda)./Lambda.*exp(-L));
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
            if inputParam.Normalization == 1
                xi = (x1'*gg)/(gg'*gg);
            else
                xi = 1;
            end
            NMSE (k,t,tt) = norm(x1 - xi*gg,2)^2/norm(gg,2)^2;
    
            sparsity(k,t,tt) = mean(piL); 
            NMSEz(k,t,tt) = norm(y - A*x1,2)^2/norm(A*gg,2)^2;
            state.eta1(k,t) = eta1;
            state.gamma1(k,t,tt) = gamma1;
            state.gamma2(k,t,tt) = gamma2;
            if any(isnan(x1))
                error('NaN occurs!')
            end
            if norm(x1-xold,'fro')^2/norm(xold,'fro')^2 < 1e-8 && convBreaker
                break;
            end
        end
     %% M step
        
        %%%% EM Maximization method 1
    %     d = ( gammaOmega*s.^2./(gammaOmega*s.^2 + gamma2*ones(rA,1)) );
    
        gammaOmega = 1./(mean((y - A*x2).*conj(y - A*x2)) + muz);
    %     Lambda = max(mean(piL),0); 
    %     gammaPrior = 1./( ( (muPrior - m).*conj(muPrior - m) + 1./gammaV ) );
    %     muPrior = m;
        
        %%%% EM Maximization
        Lambda = max(mean(piL),0);
        gammaPrior = Lambda./mean(piL.*( (muPrior - m).*conj(muPrior - m) + 1./gammaV));
        muPrior = mean(piL.*m)/Lambda;
    %     muPrior = m;
        if norm(x1 - xoldOut,'fro')/norm(xoldOut,'fro') < 1e-8 && convBreaker
            break;
        end
        %% record EM's M-step result
    %     state.mu(t) = norm(muPrior,'fro');
        state.gammaP(t,tt) = mean(Lambda./gammaPrior);
        state.lambda(t,tt) = mean(Lambda);
        state.gammaOmega(t,tt) = gammaOmega;
        state.muP(:,t,tt) = muPrior;
        fprintf('The %d-th inner iteration of the %d-th outer iteration is over.\n',t,tt)
    end
        if tt == EMOutIter
            break;
        end
        xUAVOld = xUAV; Aold = A; TdelayOld = mean(Tdelay);
        gradXUAV = -2*real([(y-A*x1)'*(Hdx*x1);(y-A*x1)'*(Hdy*x1);(y-A*x1)'*(Hdz*x1)]);
        dX = -gradXUAV/norm(gradXUAV,2);
        for nLS = 1:20
            xUAV = xUAVOld + armijoBeta^nLS*dX;
        %             Tdelay = TdelayOld + armijoBeta^nLS*dTdelay;
            [A,sMDeriv] = sensingMatrixGenWPilot(xUAV,xR,subcarrInd,pilot,Tdelay);
            nrm = norm(A,2);
            A = A/nrm;
            A = A(:,1:N);
            if norm(y-A*x2,2)^2 <= ...
                    norm(y-Aold*x2,2)^2 + armijoSigma*armijoBeta^nLS*(gradXUAV'*dX)
                break;
            end
        end
    
        %     Hdx = sMDeriv.Phidx/nrm;
        %     Hdy = sMDeriv.Phidy/nrm;
        %     Hdz = sMDeriv.Phidz/nrm;
        Hddelay = sMDeriv.PhiddelaywLoS/nrm;
        Hddelay = Hddelay(:,1:N,:);
        %     xUAV = xUAV - stepSize*gradxUAV;
        
        nSubSize = M/Nr;
        gradTdelay = zeros(Nr,1);
        for nrr = 1:Nr
            gradTdelay(nrr) = -2*real((y((nrr-1)*nSubSize+1 : nrr*nSubSize,:) - A((nrr-1)*nSubSize+1 : nrr*nSubSize,:)*x2)'*(Hddelay(:,:,nrr)*x2));
        %         hessianTdelay = 2*vH*sum(real(Hdddd(:,:,nrr).*conj(Hbar((nrr-1)*nSubSize+1 : nrr*nSubSize,:) - HPost((nrr-1)*nSubSize+1 : nrr*nSubSize,:))),'all') + ...
        %             2*real(Hdddd2(nrr));
        %             Tdelay(nrr) = Tdelay(nrr) -...
        %                 stepSize*hessianTdelay^(-1)*gradTdelay;
        %         Tdelay(nrr) = Tdelay(nrr) - stepSize*gradTdelay(nrr);
        end
        
        ddelay = -sum(gradTdelay)/(2*pi*fc/c0)^2;
        for nLS = 1:20
            Tdelay = TdelayOld + armijoBeta^nLS*ddelay;
        %             Tdelay = TdelayOld + armijoBeta^nLS*dTdelay;
            Tdelay = Tdelay*ones(Nr,1);
            [A,sMDeriv] = sensingMatrixGenWPilot(xUAV,xR,subcarrInd,pilot,Tdelay);
            nrm = norm(A,2);
            A = A/nrm;
            A = A(:,1:N);
            if norm(y-A*x2,2)^2 <= ...
                    norm(y-Aold*x2,2)^2 + armijoSigma*armijoBeta^nLS*(gradTdelay'*ddelay)
                break;
            end
        end
    
    Hdx = sMDeriv.Phidx/nrm;
    Hdy = sMDeriv.Phidy/nrm;
    Hdz = sMDeriv.Phidz/nrm;

    Hdx = Hdx(:,1:N);
    Hdy = Hdy(:,1:N);
    Hdz = Hdz(:,1:N);
    
    state.xUAVRec(:,tt)   = xUAV;
    state.TdelayRec(:,tt) = Tdelay;
end
state.NMSE = NMSE;
state.NMSEz = NMSEz;
state.sparsity = sparsity;
end