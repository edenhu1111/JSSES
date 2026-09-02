function [xest,dPost, state] = MVVBIOG(yp,yd,Phip,Phid, inputParam)
%%   Multi-view variation bayesian inference
%   y: measurements, Phi: measurement matrix, sigma: noise std. The
%   algorithm is executed till ||y - Phi x|| < eps if Nmax = 0 else Nmax iterations. Returns sparse x and
%   time elapsed
% tic
global gg dd
gammaDF = @(x,a,b)(x.^(a-1).*exp(-b.*x).*b^(a)./gamma(a));
%% Initialization
constell = inputParam.constell;
a    = 1e0;  b = 1e0;
abar = 1e0;   bbar = 1e-6;
c = 1e-3;   d = 1e-3;

Mp = size(Phip,1);
Md = size(Phid,1);
N = size(Phip,2);
T = size(Phip,3);
eps = 1e-6;
convBreak = inputParam.convBreak;
% gammaOmegaP = inputParam.gammaOmega*ones(T,1);
% gammaOmegaD = gammaOmegaP*(Md/Mp);
gammaOmega = inputParam.gammaOmega*ones(T,1);
gammaOmegaD = gammaOmega;
gammaOmegaP = gammaOmega;
Niter = inputParam.Niter;
NiterA = inputParam.NiterA;
dataNum = inputParam.dataNum;
JESD = inputParam.JESD;

xUAV = inputParam.xUAV;
% xR = inputParam.xR;
dataSubcarrInd = inputParam.dataSubcarrInd;
pilotSubcarrInd = inputParam.pilotSubcarrInd;
sMGParam = inputParam.sMGParam;
xGrid = sMGParam.grid.x;
yGrid = sMGParam.grid.y;
deltaX = zeros(size(xGrid));
deltaY = zeros(size(yGrid));
xxGrid = xGrid;
yyGrid = yGrid;


maxStep = (yGrid(2) - yGrid(1))/2;
derivP = sMGParam.derivP;
derivD = sMGParam.derivD;
PhipSearch = Phip;
derivPSearch = derivP;
PhidSearch = Phid;
derivDSearch = derivD;

lambda0 = inputParam.lambda0;
lambdaS = inputParam.lambdaS;
NaNv = size(inputParam.eyeS,1);

st = lambdaS*lambda0*ones(N,T);
% st = zeros(N,T);

st(1,:) = ones(1,T);
stbar = st;
% s = ones(N,1);
lambdaIn = lambdaS*lambda0*ones(N-1,T);

rho = (a.*st + abar.*(1-st))./(b.*st + bbar.*(1-st)).*ones(N,T);
% rho = ones(N,T);

if any(isnan(rho))
    rho = 100*ones(N,T);
end

xest = eps*ones(N,T);
xold = xest;
PhidPost = zeros(size(Phid));
% PhidPostSquared = zeros(size(Phid(:,:,1)'*Phid(:,:,1)));
PhidPostSquared = zeros(size(Phid,2),size(Phid,2),T);
v  = zeros(dataNum,T);
mu = zeros(dataNum,T);
dPost = zeros(dataNum,T);
% dPost = dd;
sigmaDPost = zeros(dataNum,T);
% for tt = 1:T
%     PhidPost(:,:,tt) = repmat(kron(ones(NaNv,1),dPost(:,tt)),1,N).*Phid(:,:,tt);
%     PhidPostSquared(:,:,tt) = Phid(:,:,tt)'*(repmat(kron(ones(NaNv,1),abs(dPost(:,tt)).^2 + sigmaDPost(:,tt)),1,N).*Phid(:,:,tt));
% end
lnrho = zeros(N,T);
nmse = zeros(NiterA,Niter);
gx = zeros(T,length(xGrid));
gy = gx;

% stepXOld = gxOld;
% stepYOld = gyOld;

armijoSigma = inputParam.armijoSigma;
armijoBeta = inputParam.armijoBeta;
attenFac = inputParam.attenFac;
armijoAlpha = 1;
iterArmijoMax = 1;
%% Coarse Estimation of D
if 1
    for niterA = 1:NiterA
        xold = xest;
        for tt = 1:T  
            %% update X
            Sx(:,:,tt) = inv((Phip(:,:,tt)'*Phip(:,:,tt))*gammaOmegaP(tt) +...
                 + diag(rho(:,tt)));
            xest(:,tt) = Sx(:,:,tt)*( Phip(:,:,tt)'*yp(:,tt)*gammaOmegaP(tt));
            muzP(tt) = abs(trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)'));
    
            ss = abs(diag(Sx(:,:,tt)));
            %% update rho and gammaOmega
%             aa(:,tt) = (stbar(:,tt).*a + (1-stbar(:,tt)).*abar + 1);
%             bb(:,tt) = (abs(xest(:,tt)).^2 + ss + stbar(:,tt).*b + (1-stbar(:,tt)).*bbar);
            aa(:,tt) = (st(:,tt).*a + (1-st(:,tt)).*abar + 1);
            bb(:,tt) = (abs(xest(:,tt)).^2 + ss + st(:,tt).*b + (1-st(:,tt)).*bbar);
%             rho(:,tt) =  (st(:,tt).*a + (1-st(:,tt)).*abar + 1)./...
%                 (abs(xest(:,tt)).^2 + ss + st(:,tt).*b + (1-st(:,tt)).*bbar);
            rho(:,tt) = aa(:,tt)./bb(:,tt);
    
            gammaOmegaP(tt) = (c+Mp)./(norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),2)^2 + muzP(tt) +d);
            gammaOmegaD(tt) = gammaOmegaP(tt); 
    
            lnrho(:,tt) = psi(aa(:,tt)) - log(bb(:,tt));
            %% update st
%             st(2:end,tt) = b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)).*lambdaIn(:,tt)./...
%                 (b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)).*lambdaIn(:,tt) + ...
%                 bbar^abar/gamma(abar).*exp((abar-1).*lnrho(2:end,tt) - bbar.*rho(2:end,tt)).*(1-lambdaIn(:,tt)));
            %% update D

            if (niterA == NiterA) && JESD 
                vv = abs(diag(Phid(:,:,tt)*(xest(:,tt)*xest(:,tt)' + Sx(:,:,tt))*Phid(:,:,tt)'));
                for ii = 1:dataNum
                    v(ii,tt)  = sum(vv(ii:dataNum:end));
                    mu(ii,tt) = conj(yd(ii:dataNum:end,tt)'*Phid(ii:dataNum:end,:,tt)*xest(:,tt));
                    qq = -v(ii,tt)*gammaOmegaD(tt)*abs(constell-mu(ii,tt)./v(ii,tt)).^2;
                    qq = qq - max(qq);
                    dPost(ii,tt) = sum(constell.*exp(qq))./sum(exp(qq));
                    sigmaDPost(ii,tt) = abs(sum(abs(constell).^2.*exp(qq))./...
                        sum(exp(qq)) - abs(dPost(ii,tt))^2);
                end
    %                 dPost = dd;
        %         PhidPost(:,:,tt) = sensingMatrixGenWSymbol(xUAV(:,tt),xR,dataSubcarrInd ,dPost(:,tt),0);
                PhidPost(:,:,tt) = repmat(kron(ones(NaNv,1),dPost(:,tt)),1,N).*Phid(:,:,tt);
                PhidPostSquared(:,:,tt) = Phid(:,:,tt)'*(repmat(kron(ones(NaNv,1),abs(dPost(:,tt)).^2 + sigmaDPost(:,tt)),1,N).*Phid(:,:,tt));
%                 PhidPostSquared2(:,:,tt) = Phid(:,:,tt)'*(repmat(kron(ones(NaNv,1),sigmaDPost(:,tt)),1,N).*Phid(:,:,tt));
        %         PhidPostSquared(:,:,tt) = PhidPost(:,:,tt)'*PhidPost(:,:,tt);
            end

        end
    end
%% Subgraph B
%     pist = st;
%     for tt = 1:T
% %         st(2:end,tt) = b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)).*lambdaIn(:,tt)./...
% %             (b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)).*lambdaIn(:,tt) + ...
% %             bbar^abar/gamma(abar).*exp((abar-1).*lnrho(2:end,tt) - bbar.*rho(2:end,tt)).*(1-lambdaIn(:,tt)));
%         piIn(:,tt) = (pist(2:end,tt)./lambdaIn(:,tt)) ./...
%             (pist(2:end,tt)./lambdaIn(:,tt) + (1 - pist(2:end,tt))./(1 - lambdaIn(:,tt)));
% 
%         pit2s(:,tt) = (piIn(:,tt).*lambda0 + (1-lambda0).*(1 - piIn(:,tt)))./...
%             (piIn(:,tt).*lambda0 + (1-lambda0).*(1 - piIn(:,tt)) + (1 - piIn(:,tt)));
%     end
%     pit2s(pit2s > 1-eps) = 1-eps;
%     pit2s(pit2s < eps)   = eps;
%     for tt = 1:T
%         mpInd = [1:tt-1,tt+1:T];
%         pis2t(:,tt) = (prod(pit2s(:,mpInd),2).*lambdaS)./(prod(pit2s(:,mpInd),2).*lambdaS +...
%             prod( (1-pit2s(:,mpInd)) ,2).*(1-lambdaS));
%         
%         lambdaIn(:,tt) = pis2t(:,tt).*lambda0;
%     end
% %     lambdaIn(lambdaIn > 1 - eps) = 1-eps;
% %     lambdaIn(lambdaIn < eps) = eps;
%     pis = (pis2t(:,1).*pit2s(:,1).*lambdaS)./(pis2t(:,1).*pit2s(:,1).*...
%         lambdaS + (1-pis2t(:,1)).*(1-pit2s(:,1)).*(1 - lambdaS));
%     stbar(2:end,:) = lambdaIn(:,:);
%     % Break when converging
% 
%     lambda0 = max(sum( repmat(pis,1,T).*st(2:end,:),'all')/sum(repmat(pis,1,T),'all'),eps);
% %     lambda0 = max(min(sum(lambdaS.*gammaDF(rho(2:end,:),a,b),'all')*lambda0./...
% %         (sum(lambdaS.*gammaDF(rho(2:end,:),a,b),'all')*lambda0 + sum(lambdaS.*gammaDF(rho(2:end),abar,bbar),'all').*(1 - lambda0)),...
% %         1-eps),eps);
%     lambdaS = min(max(sum(pis,'all')/numel(pis),eps),1-eps);
    
end
% st = lambdaS*lambda0*ones(N,T);
% % st = zeros(N,T);
% 
% st(1,:) = ones(1,T);
% stbar = st;
% % s = ones(N,1);
% lambdaIn = lambdaS*lambda0*ones(N-1,T);
% gammaOmega = inputParam.gammaOmega*ones(T,1);
% gammaOmegaD = gammaOmega;
% gammaOmegaP = gammaOmega;
% rho = (a.*st + abar.*(1-st))./(b.*st + bbar.*(1-st)).*ones(N,T);
% % rho = ones(N,T);
% 
% if any(isnan(rho))
%     rho = 100*ones(N,T);
% end

%% VBI Iteration (E-Step)
for niter = 1:Niter
    xoldOut = xest;
    %% Subgraph A
    for niterA = 1:NiterA
        xold = xest;
        for tt = 1:T 
            %% update D
            if (niter>1 && niterA == 1) && JESD
                vv = abs(diag(Phid(:,:,tt)*(xest(:,tt)*xest(:,tt)' + Sx(:,:,tt))*Phid(:,:,tt)'));
                for ii = 1:dataNum
                    v(ii,tt)  = sum(vv(ii:dataNum:end));
                    mu(ii,tt) = conj(yd(ii:dataNum:end,tt)'*Phid(ii:dataNum:end,:,tt)*xest(:,tt));
                    qq = -v(ii,tt)*gammaOmegaD(tt)*abs(constell-mu(ii,tt)./v(ii,tt)).^2;
                    qq = qq - max(qq);
                    dPost(ii,tt) = sum(constell.*exp(qq))./sum(exp(qq));
                    sigmaDPost(ii,tt) = abs(sum(abs(constell).^2.*exp(qq))./...
                        sum(exp(qq)) - abs(dPost(ii,tt))^2);
                end
%                 dPost = dd;
        %         PhidPost(:,:,tt) = sensingMatrixGenWSymbol(xUAV(:,tt),xR,dataSubcarrInd ,dPost(:,tt),0);
                PhidPost(:,:,tt) = repmat(kron(ones(NaNv,1),dPost(:,tt)),1,N).*Phid(:,:,tt);
        %         PhidPostSquared(:,:,tt) = PhidPost(:,:,tt)'*PhidPost(:,:,tt);
                PhidPostSquared(:,:,tt) = Phid(:,:,tt)'*(repmat(kron(ones(NaNv,1),abs(dPost(:,tt)).^2 + sigmaDPost(:,tt)),1,N).*Phid(:,:,tt));
%                 PhidPostSquared2(:,:,tt) = Phid(:,:,tt)'*(repmat(kron(ones(NaNv,1),sigmaDPost(:,tt)),1,N).*Phid(:,:,tt));
             end
            %% update X
            if JESD
                Sx(:,:,tt) = inv((Phip(:,:,tt)'*Phip(:,:,tt))*gammaOmegaP(tt) +...
                    PhidPostSquared(:,:,tt)*gammaOmegaD(tt) + diag(rho(:,tt)));


%                 SigmaTemp = ( eye(Mp + Md) + ...
%                     (repmat([gammaOmegaP(tt)*ones(Mp,1);gammaOmegaD(tt)*(kron(ones(NaNv,1),abs(dPost(:,tt)).^2 + sigmaDPost(:,tt)))],1,N).*...
%                     [Phip(:,:,tt);Phid(:,:,tt)])*(repmat(1./rho(:,tt),1,Mp + Md).* [Phip(:,:,tt)',Phid(:,:,tt)']) );
%                 Sx(:,:,tt) = ( eye(N) - (repmat(1./rho(:,tt),1,Mp + Md).* [Phip(:,:,tt)',Phid(:,:,tt)'])*(SigmaTemp\...
%                     (repmat([gammaOmegaP(tt)*ones(Mp,1);gammaOmegaD(tt)*(kron(ones(NaNv,1),abs(dPost(:,tt)).^2 + sigmaDPost(:,tt)))],1,N).*...
%                     [Phip(:,:,tt);Phid(:,:,tt)])) )*diag(1./rho(:,tt));


                

                xest(:,tt) = Sx(:,:,tt)*( Phip(:,:,tt)'*yp(:,tt)*gammaOmegaP(tt) +...
                    PhidPost(:,:,tt)'*yd(:,tt)*gammaOmegaD(tt) );
    %         if 0
%                 muzP(tt) = real(trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)')/Mp);
%                 muzD(tt) = real((xest(:,tt)'*(PhidPostSquared(:,:,tt)-PhidPost(:,:,tt)'*PhidPost(:,:,tt))*xest(:,tt)+...
%                     trace(Sx(:,:,tt)*(PhidPostSquared(:,:,tt))))/Md);
%                 muz(tt) = abs((trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)');
                    
                muzP(tt) = abs(trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)'));
                muzD(tt) = abs(xest(:,tt)'*(PhidPostSquared(:,:,tt) - PhidPost(:,:,tt)'*PhidPost(:,:,tt))*xest(:,tt)+...
                    trace(Sx(:,:,tt)*(PhidPostSquared(:,:,tt)) ));
            else
                Sx(:,:,tt) = inv((Phip(:,:,tt)'*Phip(:,:,tt))*gammaOmegaP(tt) +...
                     + diag(rho(:,tt)));
                xest(:,tt) = Sx(:,:,tt)*( Phip(:,:,tt)'*yp(:,tt)*gammaOmegaP(tt));
                muzP(tt) = abs(trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)'));
            end
            ss = abs(diag(Sx(:,:,tt)));
%             %% update rho and gammaOmega
%             aa(:,tt) = (stbar(:,tt).*a + (1-stbar(:,tt)).*abar + 1);
%             bb(:,tt) = (abs(xest(:,tt)).^2 + ss + stbar(:,tt).*b + (1-stbar(:,tt)).*bbar);
            aa(:,tt) = (st(:,tt).*a + (1-st(:,tt)).*abar + 1);
            bb(:,tt) = (abs(xest(:,tt)).^2 + ss + st(:,tt).*b + (1-st(:,tt)).*bbar);
%             rho(:,tt) =  (st(:,tt).*a + (1-st(:,tt)).*abar + 1)./...
%                 (abs(xest(:,tt)).^2 + ss + st(:,tt).*b + (1-st(:,tt)).*bbar);
            rho(:,tt) = aa(:,tt)./bb(:,tt);
            if JESD
%                 gammaOmegaP(tt) = (c+1)./(norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),2)^2/Mp + muzP(tt) +d);
%                 gammaOmegaD(tt) = (c+1)./((norm(yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt),2)^2)/Md + muzD(tt) +d);

%                 gammaOmega(tt) = (c+(Md+Mp))./((norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),2)^2+...
%                     norm(yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt),2)^2) + muz(tt) +d);
                gammaOmegaP(tt) = (c+Mp)./(norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),2)^2     + muzP(tt) + d);
                gammaOmegaD(tt) = (c+Md)./(norm(yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt),2)^2 + muzD(tt) + d);
            else
                gammaOmegaP(tt) = (c+Mp)./(norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),2)^2 + muzP(tt) +d);
                gammaOmegaD(tt) = gammaOmegaP(tt); 
            end
            lnrho(:,tt) = psi(aa(:,tt)) - log(bb(:,tt));
            %% update st
            st(2:end,tt) = b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)).*lambdaIn(:,tt)./...
                (b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)).*lambdaIn(:,tt) + ...
                bbar^abar/gamma(abar).*exp((abar-1).*lnrho(2:end,tt) - bbar.*rho(2:end,tt)).*(1-lambdaIn(:,tt)));
%             st(2:end,tt) = 1./...
%                 (1 + exp(log(1-lambdaIn(:,tt)) + (abar-1).*lnrho(2:end,tt) - bbar.*rho(2:end,tt)-...
%                 (log(lambdaIn(:,tt)) + (a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt))));
%             stbar(2:end,tt) = lambdaIn(:,tt);
    
        end
%         pist(pist>1-eps) = 1-eps;
%         pist(pist<eps) = eps;
        %% Break when converging
        % to be completed
        nmse(niterA,niter) = mean(vecnorm(xest-gg,2).^2./vecnorm(gg,2).^2);
        if norm(xest - xold,'fro')^2/norm(xold,'fro')^2 < 1e-5
            break;
        end
        %
        %%%%%%%%%%%%%%%%%
    end
    %% Subgraph B
    pist = st;
    for tt = 1:T
%         st(2:end,tt) = b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)).*lambdaIn(:,tt)./...
%             (b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)).*lambdaIn(:,tt) + ...
%             bbar^abar/gamma(abar).*exp((abar-1).*lnrho(2:end,tt) - bbar.*rho(2:end,tt)).*(1-lambdaIn(:,tt)));
        piIn(:,tt) = (pist(2:end,tt)./lambdaIn(:,tt)) ./...
            (pist(2:end,tt)./lambdaIn(:,tt) + (1 - pist(2:end,tt))./(1 - lambdaIn(:,tt)));

        pit2s(:,tt) = (piIn(:,tt).*lambda0 + (1-lambda0).*(1 - piIn(:,tt)))./...
            (piIn(:,tt).*lambda0 + (1-lambda0).*(1 - piIn(:,tt)) + (1 - piIn(:,tt)));
    end
    pit2s(pit2s > 1-eps) = 1-eps;
    pit2s(pit2s < eps)   = eps;
    for tt = 1:T
        mpInd = [1:tt-1,tt+1:T];
        pis2t(:,tt) = (prod(pit2s(:,mpInd),2).*lambdaS)./(prod(pit2s(:,mpInd),2).*lambdaS +...
            prod( (1-pit2s(:,mpInd)) ,2).*(1-lambdaS));
        
        lambdaIn(:,tt) = pis2t(:,tt).*lambda0;
    end
    lambdaIn(lambdaIn > 1 - eps) = 1-eps;
    lambdaIn(lambdaIn < eps) = eps;
    pis = (pis2t(:,1).*pit2s(:,1).*lambdaS)./(pis2t(:,1).*pit2s(:,1).*...
        lambdaS + (1-pis2t(:,1)).*(1-pit2s(:,1)).*(1 - lambdaS));
    stbar(2:end,:) = lambdaIn(:,:);
    % Break when converging

    %% M-step
    lambda0 = max(sum( repmat(pis,1,T).*st(2:end,:),'all')/sum(repmat(pis,1,T),'all'),eps);  
%     lambda0 = max(min(sum(lambdaS.*gammaDF(rho(2:end,:),a,b),'all')*lambda0./...
%         (sum(lambdaS.*gammaDF(rho(2:end,:),a,b),'all')*lambda0 + sum(lambdaS.*gammaDF(rho(2:end),abar,bbar),'all').*(1 - lambda0)),...
%         1-eps),eps);
    lambdaS = min(max(sum(pis,'all')/numel(pis),eps),1-eps);
    % Off-Grid
    indNZ = find(pis>0.5);
    if isempty(indNZ) 
        indNZ = 1:size(pis);
    end
    if (~JESD)
        % determine gradients and steps 
        for tt = 1:T
            PhipS = Phip(:,:,tt)*Sx(:,:,tt);
            PhipS = PhipS(:,2:end);
            gx(tt,:) = -2*gammaOmegaP(tt)*(real( (yp(:,tt) - Phip(:,:,tt)*xest(:,tt))'*(derivP(tt).Phidx.*repmat(xest(2:end,tt).',size(yp,1),1)) ) - ...
                real(diag( PhipS'*derivP(tt).Phidx).'));
            gy(tt,:) = -2*gammaOmegaP(tt)*(real( (yp(:,tt) - Phip(:,:,tt)*xest(:,tt))'*(derivP(tt).Phidy.*repmat(xest(2:end,tt).',size(yp,1),1)) ) - ...
                real(diag( PhipS'*derivP(tt).Phidy).'));

            fBase(tt) = gammaOmegaP(tt)*(norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),'fro')^2 + ...
                trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)'));
            fBase(tt) = real(fBase(tt));
        end
        gx = real(sum(gx,1)); gy = real(sum(gy,1));
        if 1
             alpha = min([abs(maxStep - deltaX);...
                abs(-maxStep - deltaX);...
                abs(maxStep - deltaY);...
                abs(-maxStep - deltaY)],[],'all' )./max(abs([gx,gy]),[],'all');
        end
%         stepX = -alpha*(gx - stepXOld/alpha*(gx'*gx)/(gxOld'*gxOld));
%         stepY = -alpha*(gy - stepYOld/alpha*(gy'*gy)/(gyOld'*gyOld));
        stepX = zeros(size(gx));stepY = zeros(size(gy));
        stepX(indNZ) = -alpha*gx(indNZ)*attenFac^(niter-1);
        stepY(indNZ) = -alpha*gy(indNZ)*attenFac^(niter-1);
        PhipSearchtmp = zeros(Mp,1+length(indNZ),T);
%         PhidSearchtmp = zeros(Md,1+length(indNZ),T);
        % armijo principle is adopted here the step is limited in 
        for iterArmijo = 1:iterArmijoMax
            xxGrid = xGrid + deltaX + armijoBeta^iterArmijo.*stepX;
            yyGrid = yGrid + deltaY + armijoBeta^iterArmijo.*stepY;
            sMGParam.grid.x = xxGrid(indNZ);
            sMGParam.grid.y = yyGrid(indNZ);
            for tt = 1:T
                sMGParam.xUAVinput = xUAV(:,tt);
                sMGParam.subcarrInd = pilotSubcarrInd;
                sMGParam.symbol = ones(length(pilotSubcarrInd),1);
                [PhipSearchtmp(:,:,tt),derivPSearchtmp(tt)] = sensingMatrixGenWSymbolOG(sMGParam);
                PhipSearch(:,1+indNZ,tt) = PhipSearchtmp(:,2:end,tt)*sqrt(Mp/Md)*sqrt(length(indNZ)/(N+1));
                derivPSearch(tt).Phidx(:,indNZ) = derivPSearchtmp(tt).Phidx*sqrt(Mp/Md)*sqrt(length(indNZ)/(N+1));
                derivPSearch(tt).Phidy(:,indNZ) = derivPSearchtmp(tt).Phidy*sqrt(Mp/Md)*sqrt(length(indNZ)/(N+1));

                f(tt) = gammaOmegaP(tt)*( norm(yp(:,tt) - PhipSearch(:,:,tt)*xest(:,tt),'fro')^2 + ...
                    trace(PhipSearch(:,:,tt)*Sx(:,:,tt)*PhipSearch(:,:,tt)'));
                f(tt) = real(f(tt));
            end
            if sum(f) <= sum(fBase) + armijoSigma*armijoBeta^iterArmijo*(gx*stepX' + gy*stepY') ||...
                iterArmijo == iterArmijoMax
                Phip = PhipSearch;
                derivP = derivPSearch;
%                 xGrid = xxGrid;
%                 yGrid = yyGrid;
%                 sMGParam.grid.x = xGrid;
%                 sMGParam.grid.y = yGrid;
                deltaX = deltaX + armijoBeta^iterArmijo.*stepX;
                deltaY = deltaY + armijoBeta^iterArmijo.*stepY;
                if JESD
                    for tt = 1:T
                        sMGParam.xUAVinput = xUAV(:,tt);
                        sMGParam.subcarrInd = dataSubcarrInd;
                        sMGParam.symbol = ones(length(dataSubcarrInd),1);
                        [Phid(:,:,tt),derivD(tt)] = sensingMatrixGenWSymbolOG(sMGParam);
                    end
                end
%                 gxOld = gx;
%                 gyOld = gy;
%                 stepXOld = stepX;
%                 stepYOld = stepY;
                break;
            end
        end

    else
        % determine gradients and steps 

        for tt = 1:T
            PhipS = Phip(:,:,tt)*Sx(:,:,tt);
            PhipS = PhipS(:,2:end);
            PhidS = Phid(:,:,tt)*Sx(:,:,tt);  
            PhidS = repmat(kron(ones(NaNv,1),abs(dPost(:,tt)).^2 + sigmaDPost(:,tt)),1,N-1).*PhidS(:,2:end);
%             PhidS = repmat(kron(ones(NaNv,1),abs(dPost(:,tt)).^2),1,N-1).*PhidS(:,2:end);
%             gx(tt,:) = -2*gammaOmegaP(tt)*(  real( (yp(:,tt) - Phip(:,:,tt)*xest(:,tt))'*(derivP(tt).Phidx.*repmat(xest(2:end,tt).',size(yp,1),1)))-...
%                 real( diag( PhipS'*derivP(tt).Phidx).' )   )...
%                 -2*gammaOmegaD(tt)*real( ( (yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt) )'*...
%                 (repmat(kron(ones(NaNv,1),dPost(:,tt)),1,N-1).*derivD(tt).Phidx.*repmat(xest(2:end,tt).',size(yd,1),1))) - ...
%                 ((kron(ones(NaNv,1),sigmaDPost(:,tt)).*Phid(:,:,tt)*xest(:,tt))'*(derivD(tt).Phidx.*repmat(xest(2:end,tt).',size(yd,1),1)))- ...
%                 (diag(PhidS'*derivD(tt).Phidx).')     );
            gx1(tt,:) = -2*gammaOmegaP(tt)*(  real( (yp(:,tt) - Phip(:,:,tt)*xest(:,tt))'*(derivP(tt).Phidx.*repmat(xest(2:end,tt).',size(yp,1),1)))-...
                real( diag( PhipS'*derivP(tt).Phidx).' )   );
            gx2(tt,:) = -2*gammaOmegaD(tt)*real( ( (yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt) )'*...
                (repmat(kron(ones(NaNv,1),dPost(:,tt)),1,N-1).*derivD(tt).Phidx.*repmat(xest(2:end,tt).',size(yd,1),1))) - ...
                ((kron(ones(NaNv,1),sigmaDPost(:,tt)).*Phid(:,:,tt)*xest(:,tt))'*(derivD(tt).Phidx.*repmat(xest(2:end,tt).',size(yd,1),1)))- ...
                (diag(PhidS'*derivD(tt).Phidx).')     );
            gx0(tt,:) = (gx1(tt,:) + gx2(tt,:));
%             gy(tt,:) = -2*gammaOmegaP(tt)*(  real((yp(:,tt) - Phip(:,:,tt)*xest(:,tt))'*(derivP(tt).Phidy.*repmat(xest(2:end,tt).',size(yp,1),1)))-...
%                 real( diag( PhipS'*derivP(tt).Phidy).' )  )...
%                 -2*gammaOmegaD(tt)*real( (yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt) )'*...
%                 (repmat(kron(ones(NaNv,1),dPost(:,tt)),1,N-1).*derivD(tt).Phidy.*repmat(xest(2:end,tt).',size(yd,1),1)) - ...
%                 ((kron(ones(NaNv,1),sigmaDPost(:,tt)).*Phid(:,:,tt)*xest(:,tt))'*(derivD(tt).Phidy.*repmat(xest(2:end,tt).',size(yd,1),1)))- ...
%                 (diag(PhidS'*derivD(tt).Phidy).')     );
            gy1(tt,:) = -2*gammaOmegaP(tt)*(  real((yp(:,tt) - Phip(:,:,tt)*xest(:,tt))'*(derivP(tt).Phidy.*repmat(xest(2:end,tt).',size(yp,1),1)))-...
                real( diag( PhipS'*derivP(tt).Phidy).' )  );
            gy2(tt,:) = -2*gammaOmegaD(tt)*real( (yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt) )'*...
                (repmat(kron(ones(NaNv,1),dPost(:,tt)),1,N-1).*derivD(tt).Phidy.*repmat(xest(2:end,tt).',size(yd,1),1)) - ...
                ((kron(ones(NaNv,1),sigmaDPost(:,tt)).*Phid(:,:,tt)*xest(:,tt))'*(derivD(tt).Phidy.*repmat(xest(2:end,tt).',size(yd,1),1)))- ...
                (diag(PhidS'*derivD(tt).Phidy).')     );
            gy0(tt,:) = gy1(tt,:) + gy2(tt,:);
            fBase(tt) = gammaOmegaP(tt)*(norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),'fro')^2 +...
                trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)')) +...
                gammaOmegaD(tt)*(norm(yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt),'fro')^2 +...
                (xest(:,tt)'*(PhidPostSquared(:,:,tt) - PhidPost(:,:,tt)'*PhidPost(:,:,tt))*xest(:,tt)+...
                    trace(Sx(:,:,tt)*(PhidPostSquared(:,:,tt)) )) );

            fBase(tt) = real(fBase(tt));
        end
          gx = real(sum(gx0,1)); gy = real(sum(gy0,1));
        if max([gx,gy],[],'all') > maxStep
%         if 1
%             alpha = min(abs(maxStep - deltaX),[],'all')./max(abs([gx,gy]),[],'all');
            alpha = min([abs(maxStep - deltaX);...
                abs(-maxStep - deltaX);...
                abs(maxStep - deltaY);...
                abs(-maxStep - deltaY)],[],'all' )./max(abs([gx,gy]),[],'all');
        else
            alpha = 1;
        end
%         stepX = -alpha*(gx - stepXOld/alpha*(gx'*gx + gy'*gy)/(gxOld'*gxOld + gyOld'*gyOld));
%         stepY = -alpha*(gy - stepYOld/alpha*(gx'*gx + gy'*gy)/(gxOld'*gxOld + gyOld'*gyOld));
        stepX = zeros(size(gx));stepY = zeros(size(gy));
        stepX(indNZ) = -alpha*gx(indNZ)*attenFac^(niter-10);
        stepY(indNZ) = -alpha*gy(indNZ)*attenFac^(niter-10);
        PhipSearchtmp = zeros(Mp,1+length(indNZ),T);
        PhidSearchtmp = zeros(Md,1+length(indNZ),T);
        % armijo principle is adopted here the step is limited in 
        for iterArmijo = 1:iterArmijoMax
            xxGrid = xGrid + deltaX + armijoBeta^iterArmijo.*stepX;
            yyGrid = yGrid + deltaY + armijoBeta^iterArmijo.*stepY;
            sMGParam.grid.x = xxGrid(indNZ);
            sMGParam.grid.y = yyGrid(indNZ);
            for tt = 1:T
                sMGParam.xUAVinput = xUAV(:,tt);
                sMGParam.subcarrInd = pilotSubcarrInd;
                sMGParam.symbol = ones(length(pilotSubcarrInd),1);
                [PhipSearchtmp(:,:,tt),derivPSearchtmp(tt)] = sensingMatrixGenWSymbolOG(sMGParam);
                PhipSearch(:,1+indNZ,tt) = PhipSearchtmp(:,2:end,tt)*sqrt(Mp/Md)*sqrt(length(indNZ)/(N+1));
                derivPSearch(tt).Phidx(:,indNZ) = derivPSearchtmp(tt).Phidx*sqrt(Mp/Md)*sqrt(length(indNZ)/(N+1));
                derivPSearch(tt).Phidy(:,indNZ) = derivPSearchtmp(tt).Phidy*sqrt(Mp/Md)*sqrt(length(indNZ)/(N+1));

                sMGParam.subcarrInd = dataSubcarrInd;
                sMGParam.symbol = ones(length(dataSubcarrInd),1);
                [PhidSearchtmp(:,:,tt),derivDSearchtmp(tt)] = sensingMatrixGenWSymbolOG(sMGParam);

                PhidSearch(:,1+indNZ,tt) = PhidSearchtmp(:,2:end,tt)*sqrt(length(indNZ)/(N+1));
                derivDSearch(tt).Phidx(:,indNZ) = derivDSearchtmp(tt).Phidx*sqrt(length(indNZ)/(N+1));
                derivDSearch(tt).Phidy(:,indNZ) = derivDSearchtmp(tt).Phidy*sqrt(length(indNZ)/(N+1));
%                 PhidPost(:,:,tt) = repmat(kron(ones(NaNv,1),dPost(:,tt)),1,N).*Phid(:,:,tt);
        %         PhidPostSquared(:,:,tt) = PhidPost(:,:,tt)'*PhidPost(:,:,tt);
%                 PhidPostSquared(:,:,tt) = Phid(:,:,tt)'*(repmat(kron(ones(NaNv,1),abs(dPost(:,tt)).^2 + sigmaDPost(:,tt)),1,N).*Phid(:,:,tt));
                PhidSearchPost(:,:,tt) = repmat(kron(ones(NaNv,1),dPost(:,tt)),1,N).*PhidSearch(:,:,tt);
%                 PhidSearchPostSquare(:,:,tt) = PhidSearchPost(:,:,tt)'*PhidSearchPost(:,:,tt);
                PhidSearchPostSquare(:,:,tt) = PhidSearch(:,:,tt)'*...
                    (repmat(kron(ones(NaNv,1),abs(dPost(:,tt)).^2 + sigmaDPost(:,tt)),1,N).*PhidSearch(:,:,tt));


%                 f(tt) = gammaOmega(tt)*norm(yp(:,tt) - PhipSearch(:,:,tt)*xest(:,tt),'fro')^2 + ...
%                     gammaOmega(tt)*(  norm(yd(:,tt) - PhidSearchPost(:,:,tt)*xest(:,tt),'fro')^2 + ...
%                     real(trace(PhipSearch(:,:,tt)*Sx(:,:,tt)*PhipSearch(:,:,tt)')+...
%                     trace(Sx(:,:,tt)*(PhidSearchPostSquare(:,:,tt)) ))  );
                f(tt) = gammaOmegaP(tt)*( norm(yp(:,tt) - PhipSearch(:,:,tt)*xest(:,tt),'fro')^2 + ...
                    real(trace(PhipSearch(:,:,tt)*Sx(:,:,tt)*PhipSearch(:,:,tt)'))  )+...
                    gammaOmegaD(tt)*(  norm(yd(:,tt) - PhidSearchPost(:,:,tt)*xest(:,tt),'fro')^2 + ...
                    xest(:,tt)'*(PhidSearchPostSquare(:,:,tt)-PhidSearchPost(:,:,tt)'*PhidSearchPost(:,:,tt))*xest(:,tt)+...
                    trace(Sx(:,:,tt)*(PhidSearchPostSquare(:,:,tt)) )) ;
                f(tt) = real(f(tt));
            end
            if sum(f) <= sum(fBase) + armijoSigma*armijoAlpha*armijoBeta^iterArmijo*(gx*stepX' + gy*stepY') || ...
                iterArmijo == iterArmijoMax
                Phip     = PhipSearch;
                derivP = derivPSearch;
                Phid = PhidSearch;
                derivD = derivDSearch;
                deltaX = deltaX + armijoBeta^iterArmijo.*stepX;
                deltaY = deltaY + armijoBeta^iterArmijo.*stepY;
                break;
            end
        end
        
    end
    if JESD && ~convBreak
        dDecoded = mleDecoder(dPost,constell);
        SER(niter) = sum(dDecoded ~= dd,'all')/numel(dd);
    end
    if norm(xest - xoldOut,'fro')^2/norm(xoldOut,'fro')^2 < 1e-6 && convBreak
        break;
    end
end
state.NMSE = nmse;
state.rho = rho;
state.st = st;
state.pis = pis;
state.gammaOmegaP = gammaOmegaP;
if JESD
    state.gammaOmegaD = gammaOmegaD;
    
end
if ~convBreak && JESD
    state.SER = SER;
end
state.lambda0 = lambda0;
state.lambdaS = lambdaS;
state.xGrid = xxGrid;
state.yGrid = yyGrid;

% state.telap = toc;
end