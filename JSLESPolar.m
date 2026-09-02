function [x1,state] = JSLES(y,A,inputParam)
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
% x2: estimated x 
% state: 
% -NMSE: the NMSE of the estimated x 
% -sparsity: the average sparsity of x
%               Nr->Na->Nv->Nsub
%               2D->1D: x outer, y inner (similar to vec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aa = 1;
% bb = 1e-3;
c0 = physconst('lightspeed');

gg = inputParam.gg;

M = size(A,1);
N = size(A,2);
numBS = size(A,3);

gammaPrior = inputParam.gammaPrior;
muPrior = inputParam.muPrior;
Lambda = inputParam.Lambda*ones(N,numBS);
gammaOmega = inputParam.gammaOmega*ones(numBS,1);
niter = inputParam.niter;
EMiter = inputParam.EMiter;
convBreaker = inputParam.convBreaker;
dampFacGam = inputParam.dampFacGam;
dampFac = inputParam.dampFac;
lambda0 = inputParam.lambda0;
lambdaS = inputParam.lambdaS;

flagOG = inputParam.flagOG;
flagSL = inputParam.flagSL;
flagSV = inputParam.flagSV;

sMGParam = inputParam.sMGParam;

xUAV = inputParam.xUAV;
xR   = inputParam.xR;
ULALine = inputParam.ULALine;


% pilotSubcarrInd = sMGParam.subcarrInd;
range = sMGParam.grid.range;
azi = sMGParam.grid.azi;
deltaR = zeros(size(range));
deltaAzi = zeros(size(azi));
% xxGrid = range;
% yyGrid = azi;
deltaxUAV = zeros(size(xUAV));
deltatau = 0;



maxStepR   = (max(range) - min(range))/(sMGParam.NRange-1)/2;
maxStepazi = abs(azi(2)-azi(1))/2;

derivP = sMGParam.derivP;
ASearch = A;
derivASearch = derivP;

attenFac = 1;
iterArmijoMax = inputParam.iterArmijoMax;
armijoBeta1 = inputParam.armijoBeta1;
armijoBeta2 = inputParam.armijoBeta2;
armijoBeta3 = inputParam.armijoBeta3;

armijoSigma = inputParam.armijoSigma;
armijoRho  = inputParam.armijoRho;

NMSEAll = zeros(numBS,niter,EMiter);
% NMSEz = zeros(niter,EMiter);

sparsity = zeros(numBS,niter,EMiter);

gamma2 = mean(gammaPrior,'all')*ones(numBS,1);

gammaMax = 1e+15;
eps = 10^-10;


c = 0; d = 0;
for indBS = 1:numBS
    x1(:,indBS) = Lambda(:,indBS).*muPrior(:,indBS) + (randn(N,1)+ 1j*randn(N,1))./sqrt(gammaPrior(:,indBS)/2);
    % r1 = zeros(N,1);
    r2(:,indBS) = zeros(N,1);
    
    [UU,SS,VV] = svd(A(:,:,indBS)); 
    % rA = rank(S);
    ss=diag(SS); svdMat(indBS).rA = sum((ss/max(ss))>=1e-5);
    svdMat(indBS).U = UU(:,1:svdMat(indBS).rA);
    SS = SS(1:svdMat(indBS).rA,1:svdMat(indBS).rA);
    svdMat(indBS).V = VV(:,1:svdMat(indBS).rA);
    svdMat(indBS).s = diag(SS);
    svdMat(indBS).yTilde = diag(1./svdMat(indBS).s)*svdMat(indBS).U'*y(:,indBS);
    piL(:,indBS) = Lambda(:,indBS).*ones(N,1);
end

for t = 1:EMiter
%% E step    
    xoldOut = x1;
    for k = 1:niter
        xold = x1;
        for indBS = 1:numBS
            % Equivalent LMMSE Estimator (EXT message for MMSE estimator)
            svdMat(indBS).d = ...
                ( gammaOmega(indBS)*svdMat(indBS).s.^2./(gammaOmega(indBS)*svdMat(indBS).s.^2 + gamma2(indBS)*ones(svdMat(indBS).rA,1)) );
    
            x2(:,indBS)   = r2(:,indBS) + svdMat(indBS).V*diag(svdMat(indBS).d)*(svdMat(indBS).yTilde - svdMat(indBS).V'*r2(:,indBS));
            muz(indBS)    = sum( svdMat(indBS).s.^2 ./ (gammaOmega(indBS)*svdMat(indBS).s.^2 + gamma2(indBS)));
            eta2(indBS)   = gamma2(indBS)*N./(N-sum(svdMat(indBS).d));
            gamma1(indBS) = gamma2(indBS)*sum(svdMat(indBS).d)./(N - sum(svdMat(indBS).d));
            r1(:,indBS)   = r2(:,indBS) + ...
                (N/svdMat(indBS).rA)*svdMat(indBS).V*diag(svdMat(indBS).d/mean(svdMat(indBS).d))*(svdMat(indBS).yTilde - svdMat(indBS).V'*r2(:,indBS));
    
            gamma1(indBS) = max(min(gammaMax,gamma1(indBS)),1/gammaMax);

            % Damping
            if k > 1 && t > 1 
                r1(:,indBS) = dampFac*r1(:,indBS) + (1-dampFac)*r1Old(:,indBS);
                gamma1(indBS) = dampFacGam*(gamma1(indBS)) + (1-dampFacGam)*gamma1old(indBS);
            elseif t > 1 && k == 1
                gamma1(indBS) = dampFacGam*(gamma1(indBS)) + (1-dampFacGam)*gamma1old(indBS);
            end
            r1Old(:,indBS) = r1(:,indBS); gamma1old(indBS) = gamma1(indBS);
    
            % MMSE estimator
            gammaV(:,indBS) = gamma1(indBS) + gammaPrior(:,indBS);
            m(:,indBS) = (gamma1(indBS).*r1(:,indBS)   + ...
                muPrior(:,indBS).*gammaPrior(:,indBS))./(gamma1(indBS) + gammaPrior(:,indBS));
            % m(isinf(gammaPrior)) = 0;m(isnan(gammaPrior)) = 0;
            
            sigmaSum(:,indBS) = 1./gammaPrior(:,indBS) + 1./gamma1(indBS);
            NNum  = Lambda(:,indBS)./sigmaSum(:,indBS)/pi;
            NDeno = (1-Lambda(:,indBS)).*gamma1(indBS).*exp(abs(r1(:,indBS)).^2.*(1./sigmaSum(:,indBS) - gamma1(indBS)))/pi;
%             NNum  = Lambda(:,indBS)./sigmaSum(:,indBS).*exp(-abs(r1(:,indBS)).^2.*(1./sigmaSum(:,indBS)))/pi;
%             NDeno = (1-Lambda(:,indBS)).*gamma1(indBS).*exp(-abs(r1(:,indBS)).^2.*gamma1(indBS))/pi;
            piL(:,indBS) = NNum./(NNum + NDeno);
            x1(:,indBS) = piL(:,indBS).*m(:,indBS);
            eta1(indBS) = 1./mean( piL(:,indBS).*(abs(m(:,indBS)).^2 + 1./gammaV(:,indBS)) - abs(x1(:,indBS)).^2 );
    
            % EXT message from MMSE to LMMSE
            gamma2(indBS) = eta1(indBS) - gamma1(indBS);
            gamma2(indBS) = max(min(gammaMax,gamma2(indBS)),1/gammaMax);
            r2(:,indBS) = x1(:,indBS) + gamma1(indBS).*(x1(:,indBS) - r1(:,indBS))./(eta1(indBS) - gamma1(indBS));
            

            % Damping
            if t > 1 && k > 1
                r2(:,indBS) = dampFac*r2(:,indBS) + (1-dampFac)*r2Old(:,indBS);
                gamma2(indBS) = dampFacGam*gamma2(indBS) + (1-dampFacGam)*gamma2Old(indBS);
            elseif t > 1 && k == 1
                gamma2(indBS) = dampFacGam*gamma2(indBS) + (1-dampFacGam)*gamma2Old(indBS);
            end
            gamma2Old(indBS) = gamma2(indBS); r2Old(:,indBS) = r2(:,indBS);
    
            %% record convergence information
            if inputParam.Normalization == 1
                xi = (x1'*gg)/(gg'*gg);
            else
                xi = 1;
            end
    
            sparsity(indBS,k,t) = mean(piL(:,indBS)); 
            NMSEAll(indBS,k,t) = norm(xi*x1(:,indBS)-gg(:,indBS),2).^2./norm(gg(:,indBS),2).^2;
%             NMSEz(k,t) = norm(y - A*x1,2)^2/norm(A*gg,2)^2;
            state.eta1  (indBS,k,t)   = eta1(indBS);
            state.gamma1(indBS,k,t) = gamma1(indBS);
            state.gamma2(indBS,k,t) = gamma2(indBS);

        end
        % EM update (Inner)
        for indBS = 1:numBS
            gammaOmega(indBS) = ...
                (c + M)./(norm(y(:,indBS) - A(:,:,indBS)*x2(:,indBS),'fro')^2 + muz(indBS) + d);
        end
        gammaPrior = mean(piL,2)./mean( piL.*( ( (muPrior - m).*conj(muPrior - m) + 1./gammaV ) ),2)*ones(1,numBS);

%         gammaPrior = mean(piL,'all')./mean( piL.*( ( (muPrior - m).*conj(muPrior - m) + 1./gammaV ) ),'all')*ones(N,numBS);
        % Break when x1 converges
        if norm(x1-xold,'fro')/norm(xold,'fro') < 1e-3
            break;
        end
         %% Fusion in MMSE module
        if ~flagSV && t > 1
            for indBS = 1:numBS
                
                num = (gammaPrior(:,indBS) + gamma1(indBS))./gammaPrior(:,indBS).*...
                    exp(-gamma1(indBS).^2./(gammaPrior(:,indBS) + gamma1(indBS)).*abs(r1(:,indBS)).^2);
            %         deno = 1;
                piIn(:,indBS) = 1./(1 + num);
            
                pit2s(:,indBS) = (piIn(:,indBS).*lambda0 + (1-lambda0).*(1 - piIn(:,indBS)))./...
                    (piIn(:,indBS).*lambda0 + (1-lambda0).*(1 - piIn(:,indBS)) + (1 - piIn(:,indBS)));
            end
            pit2s(pit2s > 1-eps) = 1-eps;
            pit2s(pit2s < eps)   = eps;
            for indBS = 1:numBS
                mpInd = [1:indBS-1,indBS+1:numBS];
                pis2t(:,indBS) = (prod(pit2s(:,mpInd),2).*lambdaS)./(prod(pit2s(:,mpInd),2).*lambdaS +...
                    prod( (1-pit2s(:,mpInd)) ,2).*(1-lambdaS));
                
                Lambda(:,indBS) = pis2t(:,indBS).*lambda0;
            end
            Lambda(Lambda > 1 - eps) = 1-eps;
            Lambda(Lambda < eps) = eps;
            pis = (pis2t(:,1).*pit2s(:,1).*lambdaS)./(pis2t(:,1).*pit2s(:,1).*...
                lambdaS + (1-pis2t(:,1)).*(1-pit2s(:,1)).*(1 - lambdaS));
            st = (Lambda.*piIn)./(Lambda.*piIn + (1-Lambda).*(1-piIn));
        end
    
        if ~flagSV && t > 1
%             lambda0 = max(sum( repmat(pis,1,numBS).*st,'all')/sum(repmat(pis,1,numBS),'all'),5e-3);
%             lambdaS = min(max(sum(pis,'all')/numel(pis),5e-3),1-(5e-3));
            indNZ = find(pis>=0.8);
        else
            Lambda = ones(N,1)*max(mean(piL,1),0);
            indNZ = find(mean(piL,2)>0.8);
        end
    end


    NMSE(t) = mean(vecnorm(x1-gg,2).^2./vecnorm(gg,2).^2);
    %% record EM's M-step result
%     gammaPrior = mean(piL,2)./mean( piL.*( ( (muPrior - m).*conj(muPrior - m) + 1./gammaV ) ),2)*ones(1,numBS);
    for indBS = 1:numBS
        cov2(:,:,indBS) = 1/eta2(indBS);
    end

%     state.mu(t) = norm(muPrior,'fro');
    state.gammaP(:,:,t) = mean(Lambda./gammaPrior);
    state.lambda(:,:,t) = Lambda;
    state.gammaOmega(:,t) = gammaOmega;
    if any([any(isnan(x1)),(norm(x1,'fro') > 10^5)])
        error('Not converged, Iteration %d',t);
    end
    if norm(x1 - xoldOut,'fro')/norm(xoldOut,'fro') < 1e-3 && convBreaker
        break;
    end


 %% Parameter Learning
%     X = x2;
    
    if flagOG && t > 1
        for indBS = 1:numBS
            gr(indBS,:) = -2*(real( (y(:,indBS) - A(:,:,indBS)*x2(:,indBS))'*(derivP(indBS).Phidr.*repmat(x2(:,indBS).',size(y,1),1)) ) - ...
                real(diag( cov2(:,:,indBS)*A(:,:,indBS)'*derivP(indBS).Phidr).'));
            gazi(indBS,:) = -2*(real( (y(:,indBS) - A(:,:,indBS)*x2(:,indBS))'*(derivP(indBS).Phidazi.*repmat(x2(:,indBS).',size(y,1),1)) ) - ...
                real(diag( cov2(:,:,indBS)*A(:,:,indBS)'*derivP(indBS).Phidazi).'));

            fBase(indBS) = (norm(y(:,indBS) - A(:,:,indBS)*x2(:,indBS),'fro')^2 + ...
                trace(A(:,:,indBS)*cov2(:,:,indBS)*A(:,:,indBS)'));
            fBase(indBS) = real(fBase(indBS));
        end
        gr = real(sum(gr,1)); 
        gazi = real(sum(gazi,1))./(range + deltaR);

%         if t == 1
        if 1
            alphaR = maxStepR/max(abs([gr]),[],'all');
            alphaAzi = maxStepazi./max(abs([gazi]),[],'all');
            alpha = min(alphaR,alphaAzi);
        end
        stepR = zeros(size(gr));stepAzi = zeros(size(gazi));
        stepR(indNZ) = -alpha*gr(indNZ)*attenFac^(t-1);
        stepAzi(indNZ) = -alpha*gazi(indNZ)*attenFac^(t-1);

        PhipSearchtmp = zeros(M,length(indNZ),numBS);
%         ASearch = A;
        % armijo principle is adopted here the step is limited in 
        for iterArmijo = 1:iterArmijoMax
%             xxGrid = range + deltaX + armijoBeta1^iterArmijo.*stepX;
%             yyGrid = azi + deltaY + armijoBeta1^iterArmijo.*stepY;
            ddeltaR = deltaR + armijoBeta1*armijoRho^(iterArmijo-1).*stepR;
            ddeltaAzi = deltaAzi + armijoBeta1*armijoRho^(iterArmijo-1).*stepAzi;
%             sMGParam.grid.range = xxGrid(indNZ);
%             sMGParam.grid.azi = yyGrid(indNZ);
%             ddeltaX(ddeltaX > maxStep) = maxStep; ddeltaX(ddeltaX < -maxStep) = -maxStep;
%             ddeltaY(ddeltaY > maxStep) = maxStep; ddeltaY(ddeltaY < -maxStep) = -maxStep;
            rrGrid = range + ddeltaR;
            aaGrid = azi + ddeltaAzi;

%             [mergedGridOut] = mergeGrid([rrGrid;yyGrid],x2,range,azi,0.8*maxStep);
%             rrGrid = mergedGridOut(1,:);
%             yyGrid = mergedGridOut(2,:);
            sMGParam.grid.range = rrGrid(indNZ);
            sMGParam.grid.azi = aaGrid(indNZ);
%             if flagSL
%                 xUAVInput = xUAV + deltaxUAV + armijoBeta2*armijoRho^(iterArmijo-1).*[stepxU;stepyU;stepzU];
%                 tauInput = deltatau + armijoBeta3*armijoRho^(iterArmijo-1).*steptau;
%                 sMGParam.xUAVinput = xUAVInput;
%                 sMGParam.vecTO = tauInput;
%             end

            for indBS = 1:numBS
    %             sMGParam.xUAVinput = xUAV;
    %                 sMGParam.subcarrInd = pilotSubcarrInd;
                sMGParam.xR = xR(:,indBS);
                sMGParam.ULALine = ULALine(:,indBS);
    %                 sMGParam.symbol = ones(length(pilotSubcarrInd),1);
                [PhipSearchtmp(:,:,indBS),derivPSearchtmp(indBS)] = sensingMatrixGenWSymbolPolar(sMGParam);
                ASearch(:,indNZ,indBS) = PhipSearchtmp(:,:,indBS)*sqrt(length(indNZ));
                derivASearch(indBS).Phidr(:,indNZ) = derivPSearchtmp(indBS).Phidr*sqrt(length(indNZ));
                derivASearch(indBS).Phidazi(:,indNZ) = derivPSearchtmp(indBS).Phidazi*sqrt(length(indNZ));

                derivASearch(indBS).PhidxUAV(:,indNZ) = derivPSearchtmp(indBS).PhidxUAV*sqrt(length(indNZ));
                derivASearch(indBS).PhidyUAV(:,indNZ) = derivPSearchtmp(indBS).PhidyUAV*sqrt(length(indNZ));
                derivASearch(indBS).PhidzUAV(:,indNZ) = derivPSearchtmp(indBS).PhidzUAV*sqrt(length(indNZ));
                derivASearch(indBS).Phidtau(:,indNZ)  = derivPSearchtmp(indBS).Phidtau *sqrt(length(indNZ));


                f(indBS) = ( norm(y(:,indBS) - ASearch(:,:,indBS)*x2(:,indBS),'fro')^2 + ...
                    trace(ASearch(:,:,indBS)*cov2(:,:,indBS)*ASearch(:,:,indBS)'));
                f(indBS) = real(f(indBS));
            end
            
            if sum(f) <= sum(fBase) + armijoSigma*armijoBeta1*armijoRho^(iterArmijo-1)*(gr*stepR' + gazi*stepAzi') ||...
                iterArmijo == iterArmijoMax
                A = ASearch;
                derivP = derivASearch;
%                 deltaX = zeros(size(deltaX));deltaY = zeros(size(deltaY));
                deltaR(indNZ) = rrGrid(indNZ) - range(indNZ);
                deltaAzi(indNZ) = aaGrid(indNZ) - azi(indNZ);
%                 if flagSL
%                     deltaxUAV = deltaxUAV + armijoBeta2*armijoRho^(iterArmijo-1).*[stepxU;stepyU;stepzU];
%                     deltatau  = deltatau  + armijoBeta3*armijoRho^(iterArmijo-1).*steptau;
%                 end
                if ~flagSL
                    for indBS = 1:numBS
                        
                        [UU,SS,VV] = svd(A(:,:,indBS)); 
                        % rA = rank(S);
                        ss=diag(SS); svdMat(indBS).rA = sum((ss/max(ss))>=1e-4);
                        svdMat(indBS).U = UU(:,1:svdMat(indBS).rA);
                        SS = SS(1:svdMat(indBS).rA,1:svdMat(indBS).rA);
                        svdMat(indBS).V = VV(:,1:svdMat(indBS).rA);
                        svdMat(indBS).s = diag(SS);
                        svdMat(indBS).yTilde = diag(1./svdMat(indBS).s)*svdMat(indBS).U'*y(:,indBS);
                    end
                end
                break;
            end
            ASearch = A;
            derivASearch = derivP;
        end
        
        if flagSL && t > 1
            for indBS = 1:numBS
                gxU(indBS) = -2*real(...
                    (y(:,indBS) - A(:,:,indBS)*x2(:,indBS))'*derivP(indBS).PhidxUAV*x2(:,indBS) -...
                    trace(derivP(indBS).PhidxUAV*cov2(:,:,indBS)*A(:,:,indBS)'));
                gyU(indBS) = -2*real(...
                    (y(:,indBS) - A(:,:,indBS)*x2(:,indBS))'*derivP(indBS).PhidyUAV*x2(:,indBS) -...
                    trace(derivP(indBS).PhidyUAV*cov2(:,:,indBS)*A(:,:,indBS)'));
                gzU(indBS) = -2*real(...
                    (y(:,indBS) - A(:,:,indBS)*x2(:,indBS))'*derivP(indBS).PhidzUAV*x2(:,indBS) -...
                    trace(derivP(indBS).PhidzUAV*cov2(:,:,indBS)*A(:,:,indBS)'));
                gtau(indBS) = -2*real(...
                    (y(:,indBS) - A(:,:,indBS)*x2(:,indBS))'*derivP(indBS).Phidtau*x2(:,indBS) -...
                    trace(derivP(indBS).Phidtau*cov2(:,:,indBS)*A(:,:,indBS)'));
                fBase(indBS) = norm(y(:,indBS) - A(:,:,indBS)*x2(:,indBS),'fro')^2 + ...
                    trace(A(:,:,indBS)*cov2(:,:,indBS)*A(:,:,indBS)');
                fBase(indBS) = real(fBase(indBS));
            end
            gxU = sum(gxU); gyU = sum(gyU); gzU = sum(gzU); gtau = sum(gtau)/c0;
            if 1
                alphaU   = min(1./abs([gxU,gyU,gzU]),[],'all');
                alphatau = abs((5e-9)*c0./gtau);
                alpha = min(alphaU,alphatau);
            end

%             stepxU  = -alphaU*  gxU *attenFac^(t-1);
%             stepyU  = -alphaU*  gyU *attenFac^(t-1);
%             stepzU  = -alphaU*  gzU *attenFac^(t-1)*0;
%             steptau = -alphatau*gtau*attenFac^(t-1);
            stepxU  = -alpha*  gxU *1^(t-1);
            stepyU  = -alpha*  gyU *1^(t-1);
            stepzU  = -alpha*  gzU *1^(t-1)*0;
            steptau = -alpha*  gtau*1^(t-1);
%             indNZ = 1:N;
            PhipSearchtmp = zeros(M,length(indNZ),numBS);
            for iterArmijo = 1:iterArmijoMax
                xUAVInput = xUAV     + deltaxUAV + armijoBeta2*armijoRho^(iterArmijo-1).*[stepxU;stepyU;stepzU];
                tauInput  = deltatau + armijoBeta3*armijoRho^(iterArmijo-1).*steptau;
                sMGParam.xUAVinput = xUAVInput;
                sMGParam.vecTO = tauInput;
                sMGParam.grid.range = rrGrid(indNZ);
                sMGParam.grid.azi = aaGrid(indNZ);
                for indBS = 1:numBS
                    sMGParam.xR = xR(:,indBS);
                    sMGParam.ULALine = ULALine(:,indBS);
        %                 sMGParam.symbol = ones(length(pilotSubcarrInd),1);
                    [PhipSearchtmp(:,:,indBS),derivPSearchtmp(indBS)] = sensingMatrixGenWSymbolPolar(sMGParam);
                    ASearch(:,indNZ,indBS) = PhipSearchtmp(:,:,indBS)*sqrt(length(indNZ));
                    derivASearch(indBS).Phidr(:,indNZ) = derivPSearchtmp(indBS).Phidr*sqrt(length(indNZ));
                    derivASearch(indBS).Phidazi(:,indNZ) = derivPSearchtmp(indBS).Phidazi*sqrt(length(indNZ));
    
                    derivASearch(indBS).PhidxUAV(:,indNZ) = derivPSearchtmp(indBS).PhidxUAV*sqrt(length(indNZ));
                    derivASearch(indBS).PhidyUAV(:,indNZ) = derivPSearchtmp(indBS).PhidyUAV*sqrt(length(indNZ));
                    derivASearch(indBS).PhidzUAV(:,indNZ) = derivPSearchtmp(indBS).PhidzUAV*sqrt(length(indNZ));
                    derivASearch(indBS).Phidtau(:,indNZ)  = derivPSearchtmp(indBS).Phidtau *sqrt(length(indNZ));
    
    
                    f(indBS) = ( norm(y(:,indBS) - ASearch(:,:,indBS)*x2(:,indBS),'fro')^2 + ...
                        trace(ASearch(:,:,indBS)*cov2(:,:,indBS)*ASearch(:,:,indBS)'));
                    f(indBS) = real(f(indBS));
                end
      
                if sum(f) <= sum(fBase) + armijoSigma*armijoBeta2*armijoRho^(iterArmijo-1)*(gxU*stepxU + gyU*stepyU + gzU*stepzU)...
                        + armijoSigma*armijoBeta3*armijoRho^(iterArmijo-1)*(gtau*steptau) || iterArmijo == iterArmijoMax

                    A = ASearch;
                    derivP = derivASearch;
                    deltaxUAV = deltaxUAV + armijoBeta2*armijoRho^(iterArmijo-1).*[stepxU;stepyU;stepzU];
                    deltatau  = deltatau  + armijoBeta3*armijoRho^(iterArmijo-1).*steptau;

                    for indBS = 1:numBS
                        [UU,SS,VV] = svd(A(:,:,indBS)); 
                        % rA = rank(S);
                        ss=diag(SS); svdMat(indBS).rA = sum((ss/max(ss))>=1e-5);
                        svdMat(indBS).U = UU(:,1:svdMat(indBS).rA);
                        SS = SS(1:svdMat(indBS).rA,1:svdMat(indBS).rA);
                        svdMat(indBS).V = VV(:,1:svdMat(indBS).rA);
                        svdMat(indBS).s = diag(SS);
                        svdMat(indBS).yTilde = diag(1./svdMat(indBS).s)*svdMat(indBS).U'*y(:,indBS);
                    end
                    break;
                end
                ASearch = A;
                derivASearch = derivP;
            end
        end
    end

%     state.muP(:,:,t) = muPrior;
%     state.pis(:,t) = pis;
%     if norm(x1 - xOldEM,'fro').^2/norm(xOldEM,'fro').^2 < 1e-3 && convBreaker
%         break;
%     end
end
state.NMSEAll = NMSEAll;
state.NMSE = NMSE;
if ~flagSV
    state.pis = pis;
else
    state.pis = mean(piL,2);
end
state.gridRes.range = range + deltaR;
state.gridRes.azi = azi + deltaAzi;
state.gammaPrior = gammaPrior;
if flagSL
    state.deltatau = deltatau;
    state.xUAV = xUAV + deltaxUAV;
end
% state.NMSEz = NMSEz;
state.sparsity = sparsity;
end