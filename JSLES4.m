function [x1,state] = JSLES4(y,A,inputParam)
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
c0 = physconst('lightspeed');  % Speed of light

gg = inputParam.gg;

M = size(A,1);
N = size(A,2);
numBS = size(A,3);

gammaPrior = inputParam.gammaPrior*ones(N,numBS);
muPrior = inputParam.muPrior*ones(N,numBS);
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

xUAV = sMGParam.xUAVinput;
xR   = inputParam.xR;
ULALine = inputParam.ULALine;

locGT = inputParam.locGT;

% pilotSubcarrInd = sMGParam.subcarrInd;
xGrid = sMGParam.grid.x;
yGrid = sMGParam.grid.y;
deltaX = zeros(size(xGrid));
deltaY = zeros(size(yGrid));
% xxGrid = xGrid;
% yyGrid = yGrid;
deltaxUAV = zeros(size(xUAV));
deltatau = sMGParam.vecTO;

% Joint distribution
p00 = (1-lambdaS) + lambdaS*(1-lambda0)^2;
p01 = lambdaS*lambda0*(1-lambda0);
p10 = lambdaS*lambda0*(1-lambda0);
p11 = lambdaS*lambda0^2;



maxStep = (yGrid(2) - yGrid(1))/2;
derivP = sMGParam.derivP;
ASearch = A;
derivASearch = derivP;
alphamin = 0.1;
attenFac = 1;
iterArmijoMax = inputParam.iterArmijoMax;
armijoBeta1 = inputParam.armijoBeta1;
armijoBeta2 = inputParam.armijoBeta2;
armijoBeta3 = inputParam.armijoBeta3;

% alpha = alphamin;

armijoSigma = inputParam.armijoSigma;
armijoRho  = inputParam.armijoRho;

NMSEAll = zeros(numBS,niter,EMiter);
% NMSEz = zeros(niter,EMiter);

sparsity = zeros(numBS,niter,EMiter);

gamma2 = inputParam.Lambda*inputParam.gammaPrior*ones(numBS,1);

gammaMax = 1e+15;
eps = 10^-10;

% targetNum = numel(find(gg(:,1)>0));

c = 1e-10; d = 1e-10;
for indBS = 1:numBS
     x1(:,indBS) = Lambda(:,indBS).*muPrior(:,indBS);
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

        end
        % EM update (Inner)

%         gammaPrior = mean(piL,'all')./mean( piL.*( ( (muPrior - m).*conj(muPrior - m) + 1./gammaV ) ),'all')*ones(N,numBS);
        % Break when x1 converges
        if any([any(isnan(x1)),(norm(x1,'fro') > 10^3)])
            error('Not converged, Iteration %d, Inner %d',t,k);
        end

         %% Fusion in MMSE module
        if ~flagSV && t > 1
            for indBS = 1:numBS
                mpInd = [1:indBS-1,indBS+1:numBS];
                num = (gammaPrior(:,indBS) + gamma1(indBS))./gammaPrior(:,indBS).*...
                    exp(-gamma1(indBS).^2./(gammaPrior(:,indBS) + gamma1(indBS)).*abs(r1(:,indBS)).^2);
            %         deno = 1;
                piIn(:,indBS) = 1./(1 + num);
                LLR = (p11*piIn(:,indBS)+p01*(1-piIn(:,indBS)))./(p01*piIn(:,indBS) + p00*(1-piIn(:,indBS)));
                Lambda(:,mpInd) = LLR./(1+LLR);
            end
            Lambda(Lambda>1-eps) = 1-eps;
            Lambda(Lambda  <eps) = eps;
            st = (Lambda.*piIn)./(Lambda.*piIn + (1-Lambda).*(1-piIn));
            pis  = max(st,[],2);
        end

        if norm(x1-xold,'fro')/norm(xold,'fro') < 1e-3
            break;
        end
        if 0
            for indBS = 1:numBS
                gammaOmega(indBS) = ...
                    (c + M)./(norm(y(:,indBS) - A(:,:,indBS)*x2(:,indBS),'fro')^2 + muz(indBS) + d);
            end
            if ~flagSV
                gammaPrior = mean(piL,2)./mean( piL.*( ( (muPrior - m).*conj(muPrior - m) + 1./gammaV ) ),2)*ones(1,numBS);
            else
                gammaPrior = mean(piL,2)./mean( piL.*( ( (muPrior - m).*conj(muPrior - m) + 1./gammaV ) ),2)*ones(1,numBS);
%                 gammaPrior = piL./ (piL.*( ( (muPrior - m).*conj(muPrior - m) + 1./gammaV ) ));
%                 gammaPrior = ones(N,1)*mean(piL,1)./mean( piL.*( ( (muPrior - m).*conj(muPrior - m) + 1./gammaV ) ),1);
%                 gammaPrior = mean(piL,'all')./mean( piL.*( ( (muPrior - m).*conj(muPrior - m) + 1./gammaV ) ),'all')*ones(N,numBS);
                Lambda = ones(N,1)*max(mean(piL,1),0);
            end
        end

    end
    if 1
        for indBS = 1:numBS
            gammaOmega(indBS) = ...
                (c + M)./(norm(y(:,indBS) - A(:,:,indBS)*x2(:,indBS),'fro')^2 + muz(indBS) + d);
        end
        if ~flagSV
            gammaPrior = mean(piL,2)./mean( piL.*( ( (muPrior - m).*conj(muPrior - m) + 1./gammaV ) ),2)*ones(1,numBS);
        else
            gammaPrior = mean(piL,2)./mean( piL.*( ( (muPrior - m).*conj(muPrior - m) + 1./gammaV ) ),2)*ones(1,numBS);
    %                 gammaPrior = piL./ (piL.*( ( (muPrior - m).*conj(muPrior - m) + 1./gammaV ) ));
    %                 gammaPrior = ones(N,1)*mean(piL,1)./mean( piL.*( ( (muPrior - m).*conj(muPrior - m) + 1./gammaV ) ),1);
    %                 gammaPrior = mean(piL,'all')./mean( piL.*( ( (muPrior - m).*conj(muPrior - m) + 1./gammaV ) ),'all')*ones(N,numBS);
            Lambda = ones(N,1)*max(mean(piL,1),0);
        end
    end
    % EM
    if ~flagSV && t > 1
%             lambda0 = max(sum( repmat(pis,1,numBS).*st,'all')/sum(repmat(pis,1,numBS),'all'),5e-3);
%             lambdaS = min(max(sum(pis,'all')/numel(pis),5e-3),1-(5e-3));
        indNZ = find(all([pis>=0.5,mean(abs(x1).^2,2)>0.001],2));
%         indNZ = 1:N;
%             Lambda = ones(N,1)*max(mean(piL,1),0);
    else
        indNZ = find(all([max(piL,[],2)>=0.5, mean(abs(x1).^2,2)>0.001],2));
    end
%     gammaPrior = mean(piL,'all')./mean( piL.*( ( (muPrior - m).*conj(muPrior - m) + 1./gammaV ) ),'all')*ones(N,numBS);
    %% record EM's M-step result
    NMSE(t) = mean(vecnorm(x1-gg,2).^2./vecnorm(gg,2).^2);
    if ~convBreaker
        
        NMSEz(t) = 0;
        for indBS = 1:numBS
            NMSEz(t) = NMSEz(t) + norm(y(:,indBS) - A(:,:,indBS)*x2(:,indBS))^2;
        end
        NMSEz(t) = NMSEz(t)/numBS;
    end
%     [~,I] = maxk(mean(abs(x1(indNZ,:)).^2,2),targetNum);
    locEst = [xGrid(indNZ) + deltaX(indNZ);yGrid(indNZ) + deltaY(indNZ)];
    if ~convBreaker
%         RMSE(t) = sqrt(calculate_multi_target_mse(locEst, locGT));
        RMSE(t) = evaluate_offgrid_performance(locEst, mean(abs(x1(indNZ,:)).^2,2).', locGT.',  4*sqrt(2)*maxStep);
    end

%     gammaPrior = mean(piL,2)./mean( piL.*( ( (muPrior - m).*conj(muPrior - m) + 1./gammaV ) ),2)*ones(1,numBS);
    for indBS = 1:numBS
        cov2(:,:,indBS) = 1/eta2(indBS);
%         cov2(:,:,indBS) = 0;
    end

%     state.mu(t) = norm(muPrior,'fro');
%     state.gammaP(:,:,t) = mean(Lambda./gammaPrior);
%     state.lambda(:,:,t) = Lambda;
%     state.gammaOmega(:,t) = gammaOmega;

%     if norm(x1 - xoldOut,'fro')/norm(xoldOut,'fro') < 1e-4 && convBreaker
    if t > 2
        if all([dGrid <= 1e-3 , convBreaker , norm(x1 - xoldOut,'fro')/norm(xoldOut,'fro') < 1e-3])
            break;
        end
    end
    if t == EMiter
        break;
    end

 %% Parameter Learning
%     X = x2;
    
    if flagOG && t > 1
        for indBS = 1:numBS
            gx(indBS,:) = -2*(real( (y(:,indBS) - A(:,:,indBS)*x2(:,indBS))'*(derivP(indBS).Phidx.*repmat(x2(:,indBS).',size(y,1),1)) ) - ...
                real(diag( cov2(:,:,indBS)*A(:,:,indBS)'*derivP(indBS).Phidx).'))*gammaOmega(indBS);
            gy(indBS,:) = -2*(real( (y(:,indBS) - A(:,:,indBS)*x2(:,indBS))'*(derivP(indBS).Phidy.*repmat(x2(:,indBS).',size(y,1),1)) ) - ...
                real(diag( cov2(:,:,indBS)*A(:,:,indBS)'*derivP(indBS).Phidy).'))*gammaOmega(indBS);
            jx(:,:,indBS) = (derivP(indBS).Phidx.*repmat(x2(:,indBS).',size(y,1),1));
            jy(:,:,indBS) = (derivP(indBS).Phidy.*repmat(x2(:,indBS).',size(y,1),1));

            fBase(indBS) = (norm(y(:,indBS) - A(:,:,indBS)*x2(:,indBS),'fro')^2 + ...
                trace(A(:,:,indBS)*cov2(:,:,indBS)*A(:,:,indBS)'))*gammaOmega(indBS);
            fBase(indBS) = real(fBase(indBS));
        end
        gx = real(sum(gx,1)); gy = real(sum(gy,1));
        jx = sum(jx,3);       jy = sum(jy,3);
        if 1
            Dmax = maxStep;
            alphaG = Dmax./max(abs([gx,gy]),[],'all')/2;
            alpha = min([alphaG,alphamin]);
        end
        stepX = zeros(size(gx));stepY = zeros(size(gy));


        if flagSL
            for indBS = 1:numBS
                gxU(indBS) = -2*real(...
                    (y(:,indBS) - A(:,:,indBS)*x2(:,indBS))'*derivP(indBS).PhidxUAV*x2(:,indBS) -...
                    trace(derivP(indBS).PhidxUAV*cov2(:,:,indBS)*A(:,:,indBS)'))*gammaOmega(indBS);
                gyU(indBS) = -2*real(...
                    (y(:,indBS) - A(:,:,indBS)*x2(:,indBS))'*derivP(indBS).PhidyUAV*x2(:,indBS) -...
                    trace(derivP(indBS).PhidyUAV*cov2(:,:,indBS)*A(:,:,indBS)'))*gammaOmega(indBS);
                gzU(indBS) = -2*real(...
                    (y(:,indBS) - A(:,:,indBS)*x2(:,indBS))'*derivP(indBS).PhidzUAV*x2(:,indBS) -...
                    trace(derivP(indBS).PhidzUAV*cov2(:,:,indBS)*A(:,:,indBS)'))*gammaOmega(indBS);
                gtau(indBS) = -2*real(...
                    (y(:,indBS) - A(:,:,indBS)*x2(:,indBS))'*derivP(indBS).Phidtau*x2(:,indBS) -...
                    trace(derivP(indBS).Phidtau*cov2(:,:,indBS)*A(:,:,indBS)'))*gammaOmega(indBS);
%                 fBase(indBS) = norm(y(:,indBS) - A(:,:,indBS)*x2(:,indBS),'fro')^2 + ...
%                     trace(A(:,:,indBS)*cov2(:,:,indBS)*A(:,:,indBS)');
%                 fBase(indBS) = real(fBase(indBS));
            end
            gxU = sum(gxU); gyU = sum(gyU); gzU = sum(gzU); gtau = sum(gtau)/c0;
            if 1
                alphaU   = min(min(1./abs([gxU,gyU]),[],'all'),alphamin);
                alphatau = min(abs(c0*10e-9./gtau),alphamin);
                alpha = min([alphaU,alphatau,alphaG]);           
            end
            
            stepxU  = -alpha*  gxU *attenFac^(t-1);
            stepyU  = -alpha*  gyU *attenFac^(t-1);
            stepzU  = -alpha*  gzU *attenFac^(t-1)*0;
%             steptau = -sign(gtau)*30e-9*attenFac^(t-1);
            steptau = -alpha*  gtau*attenFac^(t-1)/c0;

%             stepxU  = -alphaU*  gxU *attenFac^(t-1);
%             stepyU  = -alphaU*  gyU *attenFac^(t-1);
%             stepzU  = -alphaU*  gzU *attenFac^(t-1)*0;
%             steptau = -sign(gtau)*30e-9*attenFac^(t-1);
%             steptau = -alphatau*gtau*attenFac^(t-1)/c0;
        end
        stepX(indNZ) = -alpha*gx(indNZ)*attenFac^(t-1);
        stepY(indNZ) = -alpha*gy(indNZ)*attenFac^(t-1);
%         stepX(indNZ) = -alphaG*gx(indNZ)*attenFac^(t-1);
%         stepY(indNZ) = -alphaG*gy(indNZ)*attenFac^(t-1);
        PhipSearchtmp = zeros(M,length(indNZ),numBS);
%         ASearch = A;
        % armijo principle is adopted here the step is limited in 
        for iterArmijo = 1:iterArmijoMax
%             xxGrid = xGrid + deltaX + armijoBeta1^iterArmijo.*stepX;
%             yyGrid = yGrid + deltaY + armijoBeta1^iterArmijo.*stepY;
            if flagSL
                xUAVInput = xUAV     + deltaxUAV + armijoBeta2*armijoRho^(iterArmijo-1).*[stepxU;stepyU;stepzU];
                tauInput  = deltatau + armijoBeta3*armijoRho^(iterArmijo-1).*steptau;
                sMGParam.xUAVinput = xUAVInput;
                sMGParam.vecTO = tauInput;
            end
            ddeltaX = deltaX + armijoBeta1*armijoRho^(iterArmijo-1).*stepX;
            ddeltaY = deltaY + armijoBeta1*armijoRho^(iterArmijo-1).*stepY;
%             sMGParam.grid.x = xxGrid(indNZ);
%             sMGParam.grid.y = yyGrid(indNZ);
%             ddeltaX(ddeltaX > maxStep) = maxStep; ddeltaX(ddeltaX < -maxStep) = -maxStep;
%             ddeltaY(ddeltaY > maxStep) = maxStep; ddeltaY(ddeltaY < -maxStep) = -maxStep;
            xxGrid = xGrid + ddeltaX;
            yyGrid = yGrid + ddeltaY;

%             [mergedGridOut] = mergeGrid([xxGrid;yyGrid],x2,xGrid,yGrid,maxStep);
%             xxGrid = mergedGridOut(1,:);
%             yyGrid = mergedGridOut(2,:);
            sMGParam.grid.x = xxGrid(indNZ);
            sMGParam.grid.y = yyGrid(indNZ);
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
                [PhipSearchtmp(:,:,indBS),derivPSearchtmp(indBS)] = sensingMatrixGenWSymbolOG(sMGParam);
                ASearch(:,indNZ,indBS) = PhipSearchtmp(:,:,indBS);
                derivASearch(indBS).Phidx(:,indNZ) = derivPSearchtmp(indBS).Phidx;
                derivASearch(indBS).Phidy(:,indNZ) = derivPSearchtmp(indBS).Phidy;

                derivASearch(indBS).PhidxUAV(:,indNZ) = derivPSearchtmp(indBS).PhidxUAV;
                derivASearch(indBS).PhidyUAV(:,indNZ) = derivPSearchtmp(indBS).PhidyUAV;
                derivASearch(indBS).PhidzUAV(:,indNZ) = derivPSearchtmp(indBS).PhidzUAV;
                derivASearch(indBS).Phidtau (:,indNZ)  = derivPSearchtmp(indBS).Phidtau ;


                f(indBS) = ( norm(y(:,indBS) - ASearch(:,:,indBS)*x2(:,indBS),'fro')^2 + ...
                    trace(ASearch(:,:,indBS)*cov2(:,:,indBS)*ASearch(:,:,indBS)'))*gammaOmega(indBS);
                f(indBS) = real(f(indBS));
            end
            if flagSL
                termSL = (armijoSigma*armijoBeta2*armijoRho^(iterArmijo-1)*(gxU*stepxU + gyU*stepyU + gzU*stepzU)) +...
                    armijoSigma*armijoBeta3*armijoRho^(iterArmijo-1)*(c0*gtau*steptau);
            else
                termSL = 0;
            end

            if sum(f) <= sum(fBase) + armijoSigma*armijoBeta1*armijoRho^(iterArmijo-1)*(gx*stepX' + gy*stepY') + ...
                        termSL || iterArmijo == iterArmijoMax

                ffRec(2*t-1) = sum(fBase);
                ffRec(2*t)   = sum(f);

                dGrid = sqrt(norm(armijoBeta1*armijoRho^(iterArmijo-1).*stepX,2)^2 + ...
                    norm(armijoBeta1*armijoRho^(iterArmijo-1).*stepY,2)^2)/sqrt(length(indNZ));
                A = ASearch;
                derivP = derivASearch;
%                 deltaX = zeros(size(deltaX));deltaY = zeros(size(deltaY));
                deltaX(indNZ) = xxGrid(indNZ) - xGrid(indNZ);
                deltaY(indNZ) = yyGrid(indNZ) - yGrid(indNZ);
                if flagSL
                    deltaxUAV = deltaxUAV + armijoBeta2*armijoRho^(iterArmijo-1).*[stepxU;stepyU;stepzU];
                    deltatau  = deltatau  + armijoBeta3*armijoRho^(iterArmijo-1).*steptau;
                end
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
                break;
            end
            ASearch = A;
            derivASearch = derivP;
        end
    end

%     state.muP(:,:,t) = muPrior;
%     state.pis(:,t) = pis;
%     if norm(x1 - xOldEM,'fro').^2/norm(xOldEM,'fro').^2 < 1e-3 && convBreaker
%         break;
%     end
end
state.NMSE = NMSE;
if ~convBreaker
    state.RMSE = RMSE;
    
    state.NMSEz = NMSEz;
end
if ~flagSV
    state.pis = pis;
else
    state.pis = mean(piL,2);
end
ind = find(all([state.pis>=0.5,abs(x1).^2 > 0.01*max(abs(x1).^2)],2));
state.xGrid = xGrid + deltaX;
state.yGrid = yGrid + deltaY;
state.gammaPrior = gammaPrior;
if flagSL
    state.deltatau = deltatau;
    state.xUAV = xUAV + deltaxUAV;
end
% state.NMSEz = NMSEz;
% state.loc = merge_nearby_coords([state.xGrid(indNZ);state.yGrid(indNZ)], mean(abs(x1(indNZ,:)).^2,2).', 2*sqrt(2)*maxStep);
% state.RMSEfinal = sqrt(calculate_multi_target_mse(state.loc.', locGT));

[state.RMSEfinal,state.loc] = evaluate_offgrid_performance([state.xGrid(indNZ);state.yGrid(indNZ)], mean(abs(x1(indNZ,:)).^2,2).', locGT.', 4*sqrt(2)*maxStep);
end