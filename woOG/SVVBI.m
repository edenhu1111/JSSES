function [xest,dPost, state] = SVVBI(yp,yd,Phip,Phid, inputParam)
%%   Multi-view variation bayesian inference
%   y: measurements, Phi: measurement matrix, sigma: noise std. The
%   algorithm is executed till ||y - Phi x|| < eps if Nmax = 0 else Nmax iterations. Returns sparse x and
%   time elapsed
% tic
global gg
% gammaDF = @(x,a,b)(x.^(a-1).*exp(-b.*x).*b^(a)./gamma(a));

%% Initialization
constell = inputParam.constell;
% a    = 1e0;  b = 1e1;
% abar = 1e0;   bbar = 1e-8;
SV = inputParam.SV;

a    = 1e0;  b = 1e0; 
c = 1e-6;   d = 1e-6;
abar = 1e0;   bbar = 1e-6;
% if SV
%     abar = 1e0;   bbar = 1e-6;
% else
% %     abar = 1e-50;   bbar = 1e-50;
% end


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
JESD = inputParam.JESD;
lambda0 = inputParam.lambdaS*inputParam.lambda0;

% xUAV = inputParam.xUAV;
% xR = inputParam.xR;
% dataSubcarrInd = inputParam.dataSubcarrInd;
% pilotSubcarrInd = inputParam.pilotSubcarrInd;
% lambdaIn = 0.5*ones(N-1,T);

NaNv = size(inputParam.eyeS,1);

% d = zeros(dataNum,T);
if SV
    st = lambda0*ones(N,T);
    st(1,:) = ones(1,T);
else
    st = zeros(N,T);
end

rho = (a.*st + abar.*(1-st))./(b.*st + bbar.*(1-st)).*ones(N,T);
% rho = ones(N,T);

% rho = (a)./(b).*ones(N,T);

xest = zeros(N,T);
PhidPost = zeros(size(Phid));
% PhidPostSquared = zeros(size(Phid(:,:,1)'*Phid(:,:,1)));
PhidPostSquared = zeros(size(Phid,2),size(Phid,2),T);
v  = zeros(dataNum,T);
mu = zeros(dataNum,T);
dPost = zeros(dataNum,T);
sigmaDPost = zeros(dataNum,T);
lambdaIn = lambda0*ones(N-1,T);

% lnrho = zeros(N,T);
nmse = zeros(NiterA,Niter);
%% VBI Iteration
for niter = 1:Niter
    xoldOut = xest;
    %% Subgraph A
    for niterA = 1:NiterA
        xold = xest;
        for tt = 1:T
            %% update D
            if (niter > 1 && niterA == 1) && JESD
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
        %         PhidPost(:,:,tt) = sensingMatrixGenWSymbol(xUAV(:,tt),xR,dataSubcarrInd ,dPost(:,tt),0);
                PhidPost(:,:,tt) = repmat(kron(ones(NaNv,1),dPost(:,tt)),1,N).*Phid(:,:,tt);
        %         PhidPostSquared(:,:,tt) = PhidPost(:,:,tt)'*PhidPost(:,:,tt);
                PhidPostSquared(:,:,tt) = Phid(:,:,tt)'*(repmat(kron(ones(NaNv,1),abs(dPost(:,tt)).^2 + sigmaDPost(:,tt)),1,N).*Phid(:,:,tt));
            end
            %% update X
            Sx(:,:,tt) = inv((Phip(:,:,tt)'*Phip(:,:,tt))*gammaOmegaP(tt) +...
                PhidPostSquared(:,:,tt)*gammaOmegaD(tt) + diag(rho(:,tt)));
            xest(:,tt) = Sx(:,:,tt)*( Phip(:,:,tt)'*yp(:,tt)*gammaOmegaP(tt) +...
                PhidPost(:,:,tt)'*yd(:,tt)*gammaOmegaD(tt) );
            if (niter > 1) && JESD
    %         if 0
%                 muzP(tt) = real(trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)')/Mp);
%                 muzD(tt) = real((xest(:,tt)'*(PhidPostSquared(:,:,tt)-PhidPost(:,:,tt)'*PhidPost(:,:,tt))*xest(:,tt)+...
%                     trace(Sx(:,:,tt)*(PhidPostSquared(:,:,tt))))/Md);
                muz(tt) = abs((trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)')+...
                    xest(:,tt)'*(PhidPostSquared(:,:,tt)-PhidPost(:,:,tt)'*PhidPost(:,:,tt))*xest(:,tt)+...
                    trace(Sx(:,:,tt)*(PhidPostSquared(:,:,tt)) ))/(Mp+Md));
                muzP(tt) = muz(tt);muzD(tt) = muz(tt);
            else
                muzP(tt) = abs(trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)')/Mp);
            end
            ss = abs(diag(Sx(:,:,tt)));
            %% update rho and gammaOmega
            aa(:,tt) = (st(:,tt).*a + (1-st(:,tt)).*abar + 1);
            bb(:,tt) = (abs(xest(:,tt)).^2 + ss + st(:,tt).*b + (1-st(:,tt)).*bbar);
%             rho(:,tt) =  (st(:,tt).*a + (1-st(:,tt)).*abar + 1)./...
%                 (abs(xest(:,tt)).^2 + ss + st(:,tt).*b + (1-st(:,tt)).*bbar);
            if SV
                rho(:,tt) = aa(:,tt)./bb(:,tt);
                lnrho(:,tt) = psi(aa(:,tt)) - log(bb(:,tt));
            else
                rho(:,tt) = 1./(abs(xest(:,tt)).^2 + ss);
            end
            
            if (niter > 1) && JESD
%                 gammaOmegaP(tt) = (c+1)./(norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),2)^2/Mp + muzP(tt) +d);
%                 gammaOmegaD(tt) = (c+1)./((norm(yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt),2)^2)/Md + muzD(tt) +d);

                gammaOmega(tt) = (c+1)./((norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),2)^2+...
                    norm(yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt),2)^2)/(Md + Mp) + muz(tt) +d);
                gammaOmegaP(tt) = gammaOmega(tt);gammaOmegaD(tt) = gammaOmega(tt);
            else
                gammaOmegaP(tt) = (c+1)./(norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),2)^2/Mp + muzP(tt) +d);
            end
            
            %% update st
%     
%             st(2:end,tt) = b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt))./...
%                 (b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)) + ...
%                 bbar^abar/gamma(abar).*exp((abar-1).*lnrho(2:end,tt) - bbar.*rho(2:end,tt)));
            if SV
                st(2:end,tt) = b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)).*lambdaIn(:,tt)./...
                    (b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)).*lambdaIn(:,tt) + ...
                    bbar^abar/gamma(abar).*exp((abar-1).*lnrho(2:end,tt) - bbar.*rho(2:end,tt)).*(1-lambdaIn(:,tt)));
            end
        end
        %% Break when converging
        % to be completed
        nmse(niterA,niter) = mean(vecnorm(xest-gg,2).^2./vecnorm(gg,2).^2);
        if norm(xest - xold,2)^2/norm(xold)^2 < 1e-5
            break;
        end
        %
        %%%%%%%%%%%%%%%%%
    end
    %% Subgraph B
%     for tt = 1:T
%         piIn(:,tt) = (st(2:end,tt)./lambdaIn(:,tt)) ./...
%             (st(2:end,tt)./lambdaIn(:,tt) + (1 - st(2:end,tt))./(1 - lambdaIn(:,tt)));
% 
%         pit2s(:,tt) = (1 - piIn(:,tt) - lambda0 + 2.*piIn(:,tt)*lambda0)./...
%             (2 * (1 - piIn(:,tt) - lambda0 + piIn(:,tt)*lambda0) + lambda0);
%     end
%     for tt = 1:T
%         mpInd = [1:tt-1,tt+1:T];
%         pis2t(:,tt) = prod(pit2s(:,mpInd),2)./(prod(pit2s(:,mpInd),2) + prod((1-pit2s(:,mpInd)),2));
%         lambdaIn(:,tt) = pis2t(:,tt).*lambda0;
%     end
    % Break when converging
    %% EM update
%     lambda0 = mean(st(2:end,:),'all');
    if norm(xest - xoldOut,2)^2/norm(xoldOut)^2 < 1e-5 && convBreak
        break;
    end
end
state.NMSE = nmse;
state.st = st;
state.rho = rho;
state.gammaOmega = gammaOmega;
state.lambda0 = lambda0;

% state.telap = toc;
end