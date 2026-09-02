function [freq,x1,state] = OGVAMP(y,inputParam)


gammaPrior = mean(inputParam.gammaPrior,'all');
muPrior = mean(inputParam.muPrior,'all');
Lambda = inputParam.Lambda;
gammaOmega = inputParam.gammaOmega;
niter = inputParam.niter;
EMiter = inputParam.EMiter;
convBreaker = inputParam.convBreaker;
dampFacGam = inputParam.dampFacGam;
dampFac = inputParam.dampFac;

armijoBeta = inputParam.armijoBeta;
maxArmijo = inputParam.maxArmijo;
armijoRho = inputParam.armijoRho;
armijoSigma = inputParam.armijoSigma;

% A = inputParam.A;
% A1 = inputParam.A1;
% A2 = inputParam.A2;
dim1 = inputParam.dim1; dim2 = inputParam.dim2;

if ~convBreaker
    gg = inputParam.gg;
end



% grid1 = linspace(-1,1,2*dim1+1);% [-1,1)
% grid1 = grid1(1:end-1);         
% grid2 = linspace(0,2,2*dim2+1); % [0,2)
% grid2 = grid2(1:end-1);
% grid(1,:) = kron(grid1,ones(size(grid2)));
% grid(2,:) = kron(ones(size(grid1)),grid2);
grid = inputParam.grid;

deltaGrid = zeros(size(grid));

% gg = inputParam.gg;

[A,A1,A2] = sMGen(dim1,dim2,grid);


M = size(A,1);
N = size(A,2);
% NMSE = zeros(EMiter,1);
% NMSEz = zeros(niter,EMiter);

% sparsity = zeros(niter,EMiter);

gamma2 = mean(gammaPrior);


x1 = Lambda*muPrior + (randn(N,1)+ 1j*randn(N,1))/sqrt(gammaPrior/2);
% r1 = zeros(N,1);
r2 = zeros(N,1);
eps = 10^-16;

stepGrid = zeros(2,N);

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

ASearch = A;
A1Search = A1;
A2Search = A2;

for t = 1:EMiter
%% E step    
    xoldOut = x1;
    for k = 1:niter
        %% Equivalent LMMSE Estimator(EXT message for MMSE estimator)
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


%         sparsity(k,t) = mean(piL); 
% %         NMSEz(k,t) = norm(y - A*x1,2)^2/norm(A*gg,2)^2;
%         state.eta1(k,t) = eta1;
%         state.gamma1(k,t) = gamma1;
%         state.gamma2(k,t) = gamma2;
        if norm(x1-xold,'fro')/norm(xold,'fro') < 1e-3
            break;
        end

    end
    if ~convBreaker
        NMSE(t) = norm(x1-gg)^2/norm(gg)^2;
    end
 %% M step
    gammaOmega = 1./(mean((y - A*x2).*conj(y - A*x2)) + muz);
    %%%% EM Maximization
    Lambda = max(mean(piL),0);
    gammaPrior = Lambda./mean(piL.*( (muPrior - m).*conj(muPrior - m) + 1./gammaV));
%     gammaPrior = piL./(piL.*( (muPrior - m).*conj(muPrior - m) + 1./gammaV));
    if t > 1
        f1 = -2*real((y - A*x2)'*(A1.*repmat(x2.',size(y,1),1)));
        f2 = -2*real((y - A*x2)'*(A2.*repmat(x2.',size(y,1),1)));
    %     j1 = (A1.*repmat(x2.',size(y,1),1));
    %     j2 = (A2.*repmat(x2.',size(y,1),1));
        ind = find(piL>=0.5);
    %     J = [j1(:,ind),j2(:,ind)];
        fBase = norm(y-A*x2,'fro')^2;
    
        alpha = min(1/dim1/max(abs(f1)),1/dim2/max(abs(f2)));
        alpha = min(alpha,0.2);
        stepGrid = [-alpha*f1;-alpha*f2];
    %     dd = real(J'*J + 0.05*eye(2*length(ind)))\real(J'*(y - A*x2));
    %     stepGrid(1,ind) = dd(1:length(ind)).';
    %     stepGrid(2,ind) = dd(length(ind)+1:end).';
    %     ind = find(piL<0.5);
    %     stepGrid(:,ind) = 0;
    
    
        for iterArmijo = 1:maxArmijo
            deltaGridArmijo = deltaGrid + stepGrid*armijoBeta*armijoRho^(iterArmijo-1);
            
            
            [ASearch1,A1Search1,A2Search1] = sMGen(dim1,dim2,grid(:,ind) + deltaGridArmijo(:,ind));
            ASearch(:,ind) = ASearch1;
            A1Search(:,ind) = A1Search1;
            A2Search(:,ind) = A2Search1;
    
            fnew = norm(y-ASearch*x2,'fro')^2;
        
            if fnew < fBase + ...
                    armijoSigma*(f1*stepGrid(1,:).'*armijoBeta*armijoRho^(iterArmijo-1) +...
                    f2*stepGrid(2,:).'*armijoBeta*armijoRho^(iterArmijo-1)) || iterArmijo == maxArmijo
                deltaGrid = deltaGridArmijo;
                A = ASearch; A1 = A1Search; A2 = A2Search;
                [U,S,V] = svd(A); 
                % rA = rank(S);
                ss=diag(S); rA = sum((ss/max(ss))>=1e-5);
                U = U(:,1:rA);
                S = S(1:rA,1:rA);
                V = V(:,1:rA);
                s = diag(S);
                yTilde = diag(1./s)*U'*y;
                break;
            end
        end
    end
    if inputParam.Normalization == 1
        xi = (x1'*gg)/(gg'*gg);
    else
        xi = 1;
    end
%     NMSE (t) = norm(x1 - xi*gg,2)^2/norm(gg,2)^2;
%     muPrior = m;
    if norm(x1 - xoldOut,'fro')/norm(xoldOut,'fro') < 1e-3 && convBreaker
        break;
    end
    %% record EM's M-step result
%     state.mu(t) = norm(muPrior,'fro');
%     state.gammaP(t) = mean(Lambda./gammaPrior);
%     state.lambda(t) = mean(Lambda);
%     state.gammaOmega(t) = gammaOmega;
%     state.muP(:,t) = muPrior;
end
ind = find(all([piL>=0.5,abs(x1).^2 > 0.01*max(abs(x1).^2)],2));
freq = grid(:,ind) + deltaGrid(:,ind);
if ~convBreaker
state.NMSE = NMSE;
end
state.ind = ind;
end