function [xest,dPost, state] = MVVBI(yp,yd,Phip,Phid, inputParam)
%%   Multi-view variation bayesian inference
%   y: measurements, Phi: measurement matrix, sigma: noise std. The
%   algorithm is executed till ||y - Phi x|| < eps if Nmax = 0 else Nmax iterations. Returns sparse x and
%   time elapsed
tic
global gg
%% Initialization
constell = 1/sqrt(2)*[1+1j,-1+1j,-1-1j,1-1j];
a    = 1e0;  b = 1e0;
abar = 1e0;   bbar = 1e-8;
c = 1e-8;   d = 1e-8;

Mp = size(Phip,1);
Md = size(Phid,1);
N = size(Phip,2);
T = size(Phip,3);
eps = 1e-6;
convBreak = inputParam.convBreak;
gammaOmega = inputParam.gammaOmega*ones(T,1);
Niter = inputParam.Niter;
NiterA = inputParam.NiterA;
dataNum = inputParam.dataNum;
JESD = inputParam.JESD;
% xUAV = inputParam.xUAV;
% xR = inputParam.xR;
% dataSubcarrInd = inputParam.dataSubcarrInd;
% pilotSubcarrInd = inputParam.pilotSubcarrInd;
lambda0 = inputParam.lambda0;
NaNv = size(inputParam.eyeS,1);

% d = zeros(dataNum,T);
st = lambda0*ones(N,T);
% st = zeros(N,T);

st(1,:) = ones(1,T);
% s = ones(N,1);
lambdaIn = lambda0*ones(N-1,T);

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
%% VBI Iteration (E-Step)
for niter = 1:Niter
    xoldOut = xest;
    %% Subgraph A
    for niterA = 1:NiterA
        xold = xest;
        for tt = 1:T 
            %% update D
            if (niterA > 1 || niter > 1) && JESD
                vv = abs(diag(Phid(:,:,tt)*(xest(:,tt)*xest(:,tt)' + Sx(:,:,tt))*Phid(:,:,tt)'));
                for ii = 1:dataNum
                    v(ii,tt)  = sum(vv(ii:dataNum:end));
                    mu(ii,tt) = conj(yd(ii:dataNum:end,tt)'*Phid(ii:dataNum:end,:,tt)*xest(:,tt));
        
                    dPost(ii,tt) = sum(constell.*exp(-v(ii,tt)*gammaOmega(tt)*abs(constell-mu(ii,tt)./v(ii,tt)).^2)/pi*gammaOmega(tt)*v(ii,tt))./...
                        sum(exp(-v(ii,tt)*gammaOmega(tt)*abs(constell-mu(ii,tt)./v(ii,tt)).^2)/pi*gammaOmega(tt)*v(ii,tt));
                    sigmaDPost(ii,tt) = abs(sum(abs(constell).^2.*exp(-v(ii,tt)*gammaOmega(tt)*abs(constell-mu(ii,tt)./v(ii,tt)).^2)/pi*gammaOmega(tt)*v(ii,tt))./...
                        sum(exp(-v(ii,tt)*gammaOmega(tt)*abs(constell-mu(ii,tt)./v(ii,tt)).^2)/pi*gammaOmega(tt)*v(ii,tt)) - abs(dPost(ii,tt))^2);
                end
        %         PhidPost(:,:,tt) = sensingMatrixGenWSymbol(xUAV(:,tt),xR,dataSubcarrInd ,dPost(:,tt),0);
                PhidPost(:,:,tt) = repmat(kron(ones(NaNv,1),dPost(:,tt)),1,N).*Phid(:,:,tt);
        %         PhidPostSquared(:,:,tt) = PhidPost(:,:,tt)'*PhidPost(:,:,tt);
                PhidPostSquared(:,:,tt) = Phid(:,:,tt)'*(repmat(kron(ones(NaNv,1),abs(dPost(:,tt)).^2 + sigmaDPost(:,tt)),1,N).*Phid(:,:,tt));
            end
            %% update X
            Sx(:,:,tt) = inv((Phip(:,:,tt)'*Phip(:,:,tt))*gammaOmega(tt) +...
                PhidPostSquared(:,:,tt)*gammaOmega(tt) + diag(rho(:,tt)));
            xest(:,tt) = Sx(:,:,tt)*( Phip(:,:,tt)'*yp(:,tt) + PhidPost(:,:,tt)'*yd(:,tt) )*gammaOmega(tt);
            if (niterA > 1 || niter > 1) && JESD
    %         if 0
                muz = abs((trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)') + trace(Sx(:,:,tt)*PhidPostSquared(:,:,tt)))/(Mp+Md));
            else
                muz = abs(trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)')/(Mp));
            end
            ss = abs(diag(Sx(:,:,tt)));
            %% update rho and gammaOmega
            aa(:,tt) = (st(:,tt).*a + (1-st(:,tt)).*abar + 1);
            bb(:,tt) = (abs(xest(:,tt)).^2 + ss + st(:,tt).*b + (1-st(:,tt)).*bbar);
%             rho(:,tt) =  (st(:,tt).*a + (1-st(:,tt)).*abar + 1)./...
%                 (abs(xest(:,tt)).^2 + ss + st(:,tt).*b + (1-st(:,tt)).*bbar);
            rho(:,tt) = aa(:,tt)./bb(:,tt);
            if (niterA > 1 || niter > 1) && JESD
                gammaOmega(tt) = (c+1)./((norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),2)^2 + ...
                    norm(yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt),2)^2)/(Mp+Md) + muz +d);
            else
                gammaOmega(tt) = (c+1)./(norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),2)^2/Mp + muz +d);
            end
            lnrho(:,tt) = psi(st(:,tt).*a + (1-st(:,tt)).*abar + 1) - log(abs(xest(:,tt)).^2 + ss + st(:,tt).*b + (1-st(:,tt)).*bbar);
            %% update st
    
%             st(2:end,tt) = 1./...
%                 (1 + exp(log(1-lambdaIn(:,tt)) + (abar-1).*lnrho(2:end,tt) - bbar.*rho(2:end,tt)-...
%                 (log(lambdaIn(:,tt)) + (a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt))));
            st(2:end,tt) = b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)).*lambdaIn(:,tt)./...
                (b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)).*lambdaIn(:,tt) + ...
                bbar^abar/gamma(abar).*exp(( -1).*lnrho(2:end,tt) - bbar.*rho(2:end,tt)).*(1-lambdaIn(:,tt)));
    
        end
        %% Break when converging
        % to be completed
        nmse(niterA,niter) = mean(vecnorm(xest-gg,2).^2./vecnorm(gg,2).^2);
        if norm(xest - xold,2)^2/norm(xold)^2 < 1e-5 && convBreak
            break;
        end
        %
        %%%%%%%%%%%%%%%%%
    end
    %% Subgraph B
    for tt = 1:T
        piIn(:,tt) = (st(2:end,tt)./lambdaIn(:,tt)) ./...
            (st(2:end,tt)./lambdaIn(:,tt) + (1 - st(2:end,tt))./(1 - lambdaIn(:,tt)));

        pit2s(:,tt) = (piIn(:,tt).*lambda0 + (1-lambda0).*(1 - piIn(:,tt)))./...
            (piIn(:,tt).*lambda0 + (1-lambda0).*(1 - piIn(:,tt)) + (1 - piIn(:,tt)));
    end
    for tt = 1:T
        mpInd = [1:tt-1,tt+1:T];
        pis2t(:,tt) = prod(pit2s(:,mpInd),2)./(prod(pit2s(:,mpInd),2) + prod( (1-pit2s(:,mpInd)) ,2));
        
        lambdaIn(:,tt) = pis2t(:,tt).*lambda0;
    end
    pis = (pis2t(:,1).*pit2s(:,1))./(pis2t(:,1).*pit2s(:,1) + (1-pis2t(:,1)).*(1-pit2s(:,1)));
    % Break when converging

    %% M-step
    lambda0 = max(sum( repmat(pis,1,T).*st(2:end,:),'all')/sum(repmat(pis,1,T),'all'),eps);

    if norm(xest - xoldOut,2)^2/norm(xold)^2 < 1e-5 && convBreak
        break;
    end
end
state.NMSE = nmse;
state.rho = rho;
state.st = st;
state.gammaOmega = gammaOmega;
state.telap = toc;
end