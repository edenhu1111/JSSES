function [xest,dout, state] = ALESDOG(yp,yd,Phip,Phid, inputParam)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明

global gg
gammaDF = @(x,a,b)(x.^(a-1).*exp(-b.*x).*b^(a)./gamma(a));

%% Initialization
constell = inputParam.constell;
a    = 1e0;  b = 1e0;
abar = 1e0;   bbar = 1e-6;
c = 1e-6 ;   d = 1e-6;

Mp = size(Phip,1);
Md = size(Phid,1);
N = size(Phip,2);
T = size(Phip,3);
eps = 1e-3;
convBreak = inputParam.convBreak;
gammaOmega = inputParam.gammaOmega*ones(T,1);
gammaOmegaD = gammaOmega;
gammaOmegaP = gammaOmega;
Niter = inputParam.Niter;
NiterA = inputParam.NiterA;
dataNum = inputParam.dataNum;
% JESD = 0;
% xUAV = inputParam.xUAV;
% xR = inputParam.xR;
% dataSubcarrInd = inputParam.dataSubcarrInd;
% pilotSubcarrInd = inputParam.pilotSubcarrInd;
lambda0 = inputParam.lambda0;
lambdaS = inputParam.lambdaS;
NaNv = size(inputParam.eyeS,1);

xUAV = inputParam.xUAV;
% xR = inputParam.xR;
dataSubcarrInd = inputParam.dataSubcarrInd;
pilotSubcarrInd = inputParam.pilotSubcarrInd;
sMGParam = inputParam.sMGParam;
xGrid = sMGParam.grid.x;
yGrid = sMGParam.grid.y;
maxStep = (yGrid(2) - yGrid(1))/2;
derivP = sMGParam.derivP;
derivD = sMGParam.derivD;

% d = zeros(dataNum,T);
st = lambdaS*lambda0*ones(N,T);
% st = zeros(N,T);

st(1,:) = ones(1,T);
stbar = st;
% s = ones(N,1);
lambdaIn = lambdaS*lambda0*ones(N-1,T);

rho = (a.*st + abar.*(1-st))./(b.*st + bbar.*(1-st)).*ones(N,T);
xest = zeros(N,T);
xold = xest;
PhidPost = zeros(size(Phid));
% PhidPostSquared = zeros(size(Phid(:,:,1)'*Phid(:,:,1)));
PhidPostSquared = zeros(size(Phid,2),size(Phid,2),T);
v  = zeros(dataNum,T);
mu = zeros(dataNum,T);
dPost = zeros(dataNum,T);
sigmaDPost = zeros(dataNum,T);
lnrho = zeros(N,T);
nmse = zeros(NiterA,Niter);

gx = zeros(T,length(xGrid));
gy = gx;
armijoSigma = inputParam.armijoSigma;
armijoBeta  = inputParam.armijoBeta;
armijoAlpha = 1;
iterArmijoMax = 15;
if isfield(inputParam,'decoderFlag')    
    decoderFlag = inputParam.decoderFlag;
else
    decoderFlag = 0;
end
%% VBI Iteration
for niter = 1:Niter
    xoldOut = xest;
    %% Subgraph A
    for niterA = 1:NiterA
        xold = xest;
        for tt = 1:T
            
            %% update X
            SigmaTemp = inv( eye(Mp + Md) + ...
                (repmat([gammaOmegaP(tt)*ones(Mp,1);gammaOmegaD(tt)*(kron(ones(NaNv,1),abs(dPost(:,tt)).^2 + sigmaDPost(:,tt)))],1,N).*...
                [Phip(:,:,tt);Phid(:,:,tt)])*(repmat(1./rho(:,tt),1,Mp + Md).* [Phip(:,:,tt)',Phid(:,:,tt)']) );
            Sx(:,:,tt) = ( eye(N) - (repmat(1./rho(:,tt),1,Mp + Md).* [Phip(:,:,tt)',Phid(:,:,tt)'])*SigmaTemp*...
                (repmat([gammaOmegaP(tt)*ones(Mp,1);gammaOmegaD(tt)*(kron(ones(NaNv,1),abs(dPost(:,tt)).^2 + sigmaDPost(:,tt)))],1,N).*...
                [Phip(:,:,tt);Phid(:,:,tt)]) )*diag(1./rho(:,tt));
            xest(:,tt) = Sx(:,:,tt)*( Phip(:,:,tt)'*yp(:,tt)*gammaOmegaP(tt) +...
                PhidPost(:,:,tt)'*yd(:,tt)*gammaOmegaD(tt) );
            if (niter > 1)
    %         if 0
%                 muzP(tt) = real(trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)')/Mp);
%                 muzD(tt) = real((xest(:,tt)'*(PhidPostSquared(:,:,tt)-PhidPost(:,:,tt)'*PhidPost(:,:,tt))*xest(:,tt)+...
%                     trace(Sx(:,:,tt)*(PhidPostSquared(:,:,tt))))/Md);
                muz(tt) = abs((trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)')+...
                    xest(:,tt)'*(PhidPostSquared(:,:,tt)-PhidPost(:,:,tt)'*PhidPost(:,:,tt))*xest(:,tt)+...
                    trace(Sx(:,:,tt)*(PhidPostSquared(:,:,tt)) )));
                muzP(tt) = muz(tt);muzD(tt) = muz(tt);
            else
                muzP(tt) = abs(trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)'));
            end
            ss = abs(diag(Sx(:,:,tt)));
            %% update rho and gammaOmega
            aa(:,tt) = (stbar(:,tt).*a + (1-stbar(:,tt)).*abar + 1);
            bb(:,tt) = (abs(xest(:,tt)).^2 + ss + stbar(:,tt).*b + (1-stbar(:,tt)).*bbar);
%             rho(:,tt) =  (st(:,tt).*a + (1-st(:,tt)).*abar + 1)./...
%                 (abs(xest(:,tt)).^2 + ss + st(:,tt).*b + (1-st(:,tt)).*bbar);
            rho(:,tt) = aa(:,tt)./bb(:,tt);
            if (niter > 1)
%                 gammaOmegaP(tt) = (c+1)./(norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),2)^2/Mp + muzP(tt) +d);
%                 gammaOmegaD(tt) = (c+1)./((norm(yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt),2)^2)/Md + muzD(tt) +d);

                gammaOmega(tt) = (c+(Md+Mp))./((norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),2)^2+...
                    norm(yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt),2)^2) + muz(tt) +d);
                gammaOmegaP(tt) = gammaOmega(tt);gammaOmegaD(tt) = gammaOmega(tt);
            else
                gammaOmegaP(tt) = (c+Mp)./(norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),2)^2 + muzP(tt) +d);
                gammaOmegaD(tt) = gammaOmegaP(tt); gammaOmega(tt) = gammaOmegaP(tt);
            end
            lnrho(:,tt) = psi(aa(:,tt)) - log(bb(:,tt));
            %% update st
    
%             st(2:end,tt) = 1./...
%                 (1 + exp(log(1-lambdaIn(:,tt)) + (abar-1).*lnrho(2:end,tt) - bbar.*rho(2:end,tt)-...
%                 (log(lambdaIn(:,tt)) + (a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt))));
%             st(2:end,tt) = b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)).*lambdaIn(:,tt)./...
%                 (b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)).*lambdaIn(:,tt) + ...
%                 bbar^abar/gamma(abar).*exp((abar-1).*lnrho(2:end,tt) - bbar.*rho(2:end,tt)).*(1-lambdaIn(:,tt)));
            st(2:end,tt) = b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)).*lambdaIn(:,tt)./...
                (b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)).*lambdaIn(:,tt) + ...
                bbar^abar/gamma(abar).*exp((abar-1).*lnrho(2:end,tt) - bbar.*rho(2:end,tt)).*(1-lambdaIn(:,tt)));
    
        end
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
    pis = (pis2t(:,1).*pit2s(:,1).*lambdaS)./(pis2t(:,1).*pit2s(:,1).*...
        lambdaS + (1-pis2t(:,1)).*(1-pit2s(:,1)).*(1 - lambdaS));
    stbar(2:end,:) = lambdaIn(:,:);
    %% M-step

    lambda0 = max(min(sum(lambdaS.*gammaDF(rho(2:end,:),a,b),'all')*lambda0./...
        (sum(lambdaS.*gammaDF(rho(2:end,:),a,b),'all')*lambda0 + sum(lambdaS.*gammaDF(rho(2:end),abar,bbar),'all').*(1 - lambda0)),...
        1-eps),eps);
    lambdaS = min(max(sum(pis,'all')/numel(pis),eps),1-eps);
    if niter == 1
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
        if niter == 1
            alpha = maxStep./max([gx,gy],[],'all');
        end
        stepX = -alpha*(gx);
        stepY = -alpha*(gy);
        
        % armijo principle is adopted here the step is limited in 
        for iterArmijo = 1:iterArmijoMax
            xxGrid = xGrid + armijoBeta^iterArmijo.*stepX;
            yyGrid = yGrid + armijoBeta^iterArmijo.*stepY;
            sMGParam.grid.x = xxGrid;
            sMGParam.grid.y = yyGrid;
            for tt = 1:T
                sMGParam.xUAVinput = xUAV(:,tt);
                sMGParam.subcarrInd = pilotSubcarrInd;
                sMGParam.symbol = ones(length(pilotSubcarrInd),1);
                [PhipSearch(:,:,tt),derivPSearch(tt)] = sensingMatrixGenWSymbolOG(sMGParam);
                PhipSearch(:,:,tt) = PhipSearch(:,:,tt)*sqrt(Mp/Md);
                derivPSearch(tt).Phidx = derivPSearch(tt).Phidx*sqrt(Mp/Md);
                derivPSearch(tt).Phidy = derivPSearch(tt).Phidy*sqrt(Mp/Md);

                f(tt) = gammaOmegaP(tt)*( norm(yp(:,tt) - PhipSearch(:,:,tt)*xest(:,tt),'fro')^2 + ...
                    trace(PhipSearch(:,:,tt)*Sx(:,:,tt)*PhipSearch(:,:,tt)'));
                f(tt) = real(f(tt));
            end
            if sum(f) <= sum(fBase) + armijoSigma*armijoBeta^iterArmijo*(gx*stepX' + gy*stepY')
                xGrid = xxGrid;
                yGrid = yyGrid;
                Phip = PhipSearch;
                derivP = derivPSearch;
                sMGParam.grid.x = xGrid;
                sMGParam.grid.y = yGrid;
                for tt = 1:T
                    sMGParam.xUAVinput = xUAV(:,tt);
                    sMGParam.subcarrInd = dataSubcarrInd;
                    sMGParam.symbol = ones(length(dataSubcarrInd),1);
                    [Phid(:,:,tt),derivD(tt)] = sensingMatrixGenWSymbolOG(sMGParam);
                end
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
            gx(tt,:) = -2*gammaOmegaP(tt)*(  real((yp(:,tt) - Phip(:,:,tt)*xest(:,tt))'*(derivP(tt).Phidx.*repmat(xest(2:end,tt).',size(yp,1),1)))-...
                real( diag( PhipS'*derivP(tt).Phidx).' )   )...
                -2*gammaOmegaD(tt)*real( ( (yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt) )'*...
                (repmat(kron(ones(NaNv,1),dPost(:,tt)),1,N-1).*derivD(tt).Phidx.*repmat(xest(2:end,tt).',size(yd,1),1))) - ...
                ((kron(ones(NaNv,1),sigmaDPost(:,tt)).*Phid(:,:,tt)*xest(:,tt))'*(derivD(tt).Phidx.*repmat(xest(2:end,tt).',size(yd,1),1)))- ...
                (diag(PhidS'*derivD(tt).Phidx).')     );

            gy(tt,:) = -2*gammaOmegaP(tt)*(  real((yp(:,tt) - Phip(:,:,tt)*xest(:,tt))'*(derivP(tt).Phidy.*repmat(xest(2:end,tt).',size(yp,1),1)))-...
                real( diag( PhipS'*derivP(tt).Phidy).' )  )...
                -2*gammaOmegaD(tt)*real( (yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt) )'*...
                (repmat(kron(ones(NaNv,1),dPost(:,tt)),1,N-1).*derivD(tt).Phidy.*repmat(xest(2:end,tt).',size(yd,1),1)) - ...
                ((kron(ones(NaNv,1),sigmaDPost(:,tt)).*Phid(:,:,tt)*xest(:,tt))'*(derivD(tt).Phidy.*repmat(xest(2:end,tt).',size(yd,1),1)))- ...
                (diag(PhidS'*derivD(tt).Phidy).')     );
            fBase(tt) = gammaOmegaP(tt)*(norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),'fro')^2 +...
                trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)')) +...
                gammaOmegaD(tt)*(norm(yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt),'fro')^2 +...
                (xest(:,tt)'*(PhidPostSquared(:,:,tt)-PhidPost(:,:,tt)'*PhidPost(:,:,tt))*xest(:,tt)+...
                    trace(Sx(:,:,tt)*(PhidPostSquared(:,:,tt)) )) );

%             gx(tt,:) = -2*gammaOmega(tt)*(  real((yp(:,tt) - Phip(:,:,tt)*xest(:,tt))'*(derivP(tt).Phidx.*repmat(xest(2:end,tt).',size(yp,1),1)))-...
%                 real( diag( PhipS'*derivP(tt).Phidx).' )  )...
%                 -2*gammaOmega(tt)*real( ( (yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt) )'*...
%                 (repmat(kron(ones(NaNv,1),dPost(:,tt)),1,N-1).*derivD(tt).Phidx.*repmat(xest(2:end,tt).',size(yd,1),1))) - ...
%                 (diag(PhidS'*derivD(tt).Phidx).')     );
% 
%             gy(tt,:) = -2*gammaOmega(tt)*(  real((yp(:,tt) - Phip(:,:,tt)*xest(:,tt))'*(derivP(tt).Phidy.*repmat(xest(2:end,tt).',size(yp,1),1)))-...
%                 real( diag( PhipS'*derivP(tt).Phidy).' )  )...
%                 -2*gammaOmega(tt)*real( (yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt) )'*...
%                 (repmat(kron(ones(NaNv,1),dPost(:,tt)),1,N-1).*derivD(tt).Phidy.*repmat(xest(2:end,tt).',size(yd,1),1)) - ...
%                 (diag(PhidS'*derivD(tt).Phidy).')     );
%             fBase(tt) = gammaOmega(tt)*norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),'fro')^2 + ...
%                 gammaOmega(tt)*(norm(yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt),'fro')^2 +...
%                 (trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)')+...
%                     trace(Sx(:,:,tt)*(PhidPostSquared(:,:,tt)) )) );

            fBase(tt) = real(fBase(tt));
        end
        gx = real(sum(gx,1)); gy = real(sum(gy,1));
        if niter == 1
            alpha = maxStep./max([gx,gy],[],'all');
        end
        stepX = -alpha*(gx);
        stepY = -alpha*(gy);
        % armijo principle is adopted here the step is limited in 
        for iterArmijo = 1:iterArmijoMax
            xxGrid = xGrid + armijoAlpha*armijoBeta^iterArmijo.*stepX;
            yyGrid = yGrid + armijoAlpha*armijoBeta^iterArmijo.*stepY;
            sMGParam.grid.x = xxGrid;
            sMGParam.grid.y = yyGrid;
            for tt = 1:T
                sMGParam.xUAVinput = xUAV(:,tt);
                sMGParam.subcarrInd = pilotSubcarrInd;
                sMGParam.symbol = ones(length(pilotSubcarrInd),1);
                [PhipSearch(:,:,tt),derivPSearch(tt)] = sensingMatrixGenWSymbolOG(sMGParam);

                sMGParam.subcarrInd = dataSubcarrInd;
                sMGParam.symbol = ones(length(dataSubcarrInd),1);
                [PhidSearch(:,:,tt),derivDSearch(tt)] = sensingMatrixGenWSymbolOG(sMGParam);
%                 PhidPost(:,:,tt) = repmat(kron(ones(NaNv,1),dPost(:,tt)),1,N).*Phid(:,:,tt);
        %         PhidPostSquared(:,:,tt) = PhidPost(:,:,tt)'*PhidPost(:,:,tt);
%                 PhidPostSquared(:,:,tt) = Phid(:,:,tt)'*(repmat(kron(ones(NaNv,1),abs(dPost(:,tt)).^2 + sigmaDPost(:,tt)),1,N).*Phid(:,:,tt));
                PhidSearchPost(:,:,tt) = repmat(kron(ones(NaNv,1),dPost(:,tt)),1,N).*PhidSearch(:,:,tt);
%                 PhidSearchPostSquare(:,:,tt) = PhidSearchPost(:,:,tt)'*PhidSearchPost(:,:,tt);
                PhidSearchPostSquare(:,:,tt) = PhidSearch(:,:,tt)'*...
                    (repmat(kron(ones(NaNv,1),abs(dPost(:,tt)).^2 + sigmaDPost(:,tt)),1,N).*PhidSearch(:,:,tt));

                PhipSearch(:,:,tt) = PhipSearch(:,:,tt)*sqrt(Mp/Md);
                derivPSearch(tt).Phidx = derivPSearch(tt).Phidx*sqrt(Mp/Md);
                derivPSearch(tt).Phidy = derivPSearch(tt).Phidy*sqrt(Mp/Md);
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
            if sum(f) <= sum(fBase) + armijoSigma*armijoAlpha*armijoBeta^iterArmijo*(gx*stepX' + gy*stepY') 
                xGrid = xxGrid;
                yGrid = yyGrid;
%                 sMGParam.grid.x = xGrid;
%                 sMGParam.grid.y = yGrid;
                Phip = PhipSearch;
                derivP = derivPSearch;
                Phid = PhidSearch;
                derivD = derivDSearch;

                break;
            end
        end
    end

    %% decoding
    for tt = 1:T
        vv = abs(diag(Phid(:,:,tt)*(xest(:,tt)*xest(:,tt)')*Phid(:,:,tt)'));
        for ii = 1:dataNum
            v(ii,tt)  = sum(vv(ii:dataNum:end));
            mu(ii,tt) = conj(yd(ii:dataNum:end,tt)'*Phid(ii:dataNum:end,:,tt)*xest(:,tt));
    
            dPost(ii,tt) = sum(constell.*exp(-v(ii,tt)*gammaOmega(tt)*abs(constell-mu(ii,tt)./v(ii,tt)).^2)/pi*gammaOmega(tt)*v(ii,tt))./...
                sum(exp(-v(ii,tt)*gammaOmega(tt)*abs(constell-mu(ii,tt)./v(ii,tt)).^2)/pi*gammaOmega(tt)*v(ii,tt));
    %             sigmaDPost(ii,tt) = abs(sum(abs(constell).^2.*exp(-v(ii,tt)*gammaOmega(tt)*abs(constell-mu(ii,tt)./v(ii,tt)).^2)/pi*gammaOmega(tt)*v(ii,tt))./...
    %                 sum(exp(-v(ii,tt)*gammaOmega(tt)*abs(constell-mu(ii,tt)./v(ii,tt)).^2)/pi*gammaOmega(tt)*v(ii,tt)) - abs(dPost(ii,tt))^2);
        end
        if decoderFlag
            dout = mleDecoder(dPost,constell);
        else
            dout = dPost;
        end
    %         PhidPost(:,:,tt) = sensingMatrixGenWSymbol(xUAV(:,tt),xR,dataSubcarrInd ,dPost(:,tt),0);
        PhidPost(:,:,tt) = repmat(kron(ones(NaNv,1),dout(:,tt)),1,N).*Phid(:,:,tt);
    %         PhidPostSquared(:,:,tt) = PhidPost(:,:,tt)'*PhidPost(:,:,tt);
        PhidPostSquared(:,:,tt) = PhidPost(:,:,tt)'*PhidPost(:,:,tt);
    end
    if norm(xest - xoldOut,'fro')^2/norm(xold,'fro')^2 < 1e-5 && convBreak
        break;
    end
end
dout = mleDecoder(dPost,constell);
state.NMSE = nmse;
state.rho = rho;
state.st = st;
state.pis = pis;
state.lambda0 = lambda0;
state.lambdaS = lambdaS;

state.gammaOmega = gammaOmega;
% state.telap = toc;
end